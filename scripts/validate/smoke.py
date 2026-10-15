#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""Drive SeisSol's end-to-end smoke tests and check what the binaries produce.

Registered with CTest through ``seissol_add_smoke_test()``. Two modes:

``run``
    Start a command and assert how it terminated, optionally matching its
    output against a regular expression. CTest cannot express "expected to
    terminate unsuccessfully" by itself: ``WILL_FAIL`` and
    ``PASS_REGULAR_EXPRESSION`` both report a failure when the child dies from
    a signal instead of returning a non-zero code. Which of the two happens
    depends on the build, since ``logError()`` maps to ``MPI_Abort()`` with MPI
    and to ``abort()`` without, and on how the MPI implementation turns an
    abort into an exit status.

``proxy``
    Run the proxy in its JSON output mode and check the result for formal
    correctness: it has to parse, carry exactly the documented field set, and
    hold finite numbers in a plausible range. Nothing here asserts anything
    about performance, so the checks hold on any machine.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import signal
import subprocess
import sys
from typing import Any

# Field set written by seissol::proxy::Aux::writeOutput() in JSON mode, and the
# constraint each value has to satisfy. Checking the set exactly (rather than
# just the fields of interest) means a field added to or dropped from the
# writer surfaces here as a deliberate decision instead of silently changing
# the format consumers see.
#
# "positive" is used where a zero would mean the measurement did not happen.
PROXY_FIELDS: dict[str, str] = {
    "name": "string",
    "cycle-source": "string",
    "time": "positive",
    "cycles": "non-negative",
    "gflop-libxsmm": "non-negative",
    "gflop-pspamm": "non-negative",
    "gflop-libxsmm-pspamm": "non-negative",
    "gflop-nz": "non-negative",
    "gflop-hw": "non-negative",
    "gib": "non-negative",
    "gib-kernel": "non-negative",
    "gflopcycle-nz": "non-negative",
    "gflopcycle-hw": "non-negative",
    "gibcycle": "non-negative",
    "gibcycle-kernel": "non-negative",
    "gflops-nz": "non-negative",
    "gflops-hw": "non-negative",
    "gibs": "non-negative",
    "gibs-kernel": "non-negative",
}

# Relations between fields that follow from how the numbers are defined, not
# from how fast the machine is. Each entry is (left, op, right, why).
PROXY_RELATIONS: list[tuple[str, str, str, str]] = [
    (
        "gflop-hw",
        ">=",
        "gflop-nz",
        "hardware flops include the padding the non-zero count leaves out",
    ),
]


# Mirrors CycleSource in src/Proxy/Cycles.h. A value outside this set means the
# two have drifted apart.
CYCLE_SOURCES = frozenset({"tsc", "cntvct", "none"})


class CheckFailed(Exception):
    """A smoke-test assertion did not hold."""


def describe_exit(returncode: int) -> str:
    """Render a return code the way a reader needs it, signals included."""
    if returncode < 0:
        try:
            name = signal.Signals(-returncode).name
        except ValueError:
            name = f"signal {-returncode}"
        return f"killed by {name}"
    return f"exit code {returncode}"


def run_command(
    command: list[str], timeout: float | None
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except FileNotFoundError as error:
        raise CheckFailed(f"cannot execute {command[0]!r}: {error}") from error
    except subprocess.TimeoutExpired as error:
        raise CheckFailed(f"timed out after {timeout}s") from error


def echo_output(proc: subprocess.CompletedProcess[str]) -> None:
    """Hand the child's output to CTest so --output-on-failure stays useful."""
    if proc.stdout:
        sys.stdout.write(proc.stdout)
    if proc.stderr:
        sys.stderr.write(proc.stderr)


# ---------------------------------------------------------------------------
# Mode: run
# ---------------------------------------------------------------------------


def check_termination(returncode: int, expect_failure: bool) -> None:
    if expect_failure:
        if returncode == 0:
            raise CheckFailed("expected an unsuccessful exit, got a successful one")
    elif returncode != 0:
        raise CheckFailed(f"terminated unsuccessfully: {describe_exit(returncode)}")


def check_output_matches(text: str, pattern: str | None) -> None:
    if pattern and re.search(pattern, text) is None:
        raise CheckFailed(f"output does not match: {pattern}")


def mode_run(args: argparse.Namespace) -> int:
    proc = run_command(args.command, args.timeout)
    echo_output(proc)
    check_termination(proc.returncode, args.expect_failure)
    check_output_matches(proc.stdout + proc.stderr, args.expect_output)
    return 0


# ---------------------------------------------------------------------------
# Mode: proxy
# ---------------------------------------------------------------------------


def extract_json_object(stdout: str) -> str:
    """Pick the JSON document out of the proxy's standard output.

    The writer emits the object as the final line. Anything the run printed
    before it is tolerated but reported, because a machine-readable mode that
    needs its output filtered is a problem in its own right.
    """
    lines = [line for line in stdout.splitlines() if line.strip()]
    if not lines:
        raise CheckFailed("no output to parse")
    candidate = lines[-1].strip()
    if not candidate.startswith("{"):
        raise CheckFailed(f"last output line is not a JSON object: {candidate!r}")
    if len(lines) > 1:
        preamble = " | ".join(lines[:-1])
        raise CheckFailed(
            "JSON mode wrote non-JSON lines to stdout before the object: " + preamble
        )
    return candidate


def _reject_constant(token: str) -> Any:
    raise CheckFailed(
        f"{token!r} is not valid JSON; a non-finite number reached the writer"
    )


def parse_proxy_json(text: str) -> dict[str, Any]:
    try:
        document = json.loads(text, parse_constant=_reject_constant)
    except json.JSONDecodeError as error:
        # Bare inf/nan is what a division by zero produces through
        # operator<<, and it is the most likely reason to land here.
        if re.search(r"[:,]\s*-?(inf|nan)\b", text):
            raise CheckFailed(
                "output contains bare 'inf'/'nan', which is not valid JSON; "
                "a non-finite number reached the writer"
            ) from error
        raise CheckFailed(f"output is not valid JSON: {error}") from error
    if not isinstance(document, dict):
        raise CheckFailed("top level of the output is not an object")
    return document


def check_proxy_fields(document: dict[str, Any], kernel: str) -> None:
    missing = sorted(set(PROXY_FIELDS) - set(document))
    unexpected = sorted(set(document) - set(PROXY_FIELDS))
    if missing:
        raise CheckFailed("missing field(s): " + ", ".join(missing))
    if unexpected:
        raise CheckFailed("unexpected field(s): " + ", ".join(unexpected))

    problems: list[str] = []
    for name, constraint in PROXY_FIELDS.items():
        value = document[name]
        if constraint == "string":
            if not isinstance(value, str):
                problems.append(f"{name}: expected a string, got {value!r}")
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            problems.append(f"{name}: expected a number, got {value!r}")
            continue
        if not math.isfinite(value):
            problems.append(f"{name}: not finite ({value})")
            continue
        if constraint == "positive" and value <= 0:
            problems.append(f"{name}: expected a positive value, got {value}")
        elif constraint == "non-negative" and value < 0:
            problems.append(f"{name}: expected a non-negative value, got {value}")

    if document.get("name") != kernel:
        problems.append(f"name: expected {kernel!r}, got {document.get('name')!r}")

    source = document.get("cycle-source")
    if source not in CYCLE_SOURCES:
        problems.append(
            f"cycle-source: {source!r} is not one of "
            + ", ".join(sorted(CYCLE_SOURCES))
        )
    elif source == "none" and document.get("cycles"):
        problems.append("cycle-source is 'none' but cycles is non-zero")

    combined = document.get("gflop-libxsmm-pspamm")
    parts = (document.get("gflop-libxsmm"), document.get("gflop-pspamm"))
    if all(isinstance(v, (int, float)) for v in (combined, *parts)):
        expected = parts[0] + parts[1]
        if not math.isclose(combined, expected, rel_tol=1e-9, abs_tol=1e-12):
            problems.append(
                f"gflop-libxsmm-pspamm ({combined}) is not the sum of "
                f"gflop-libxsmm and gflop-pspamm ({expected})"
            )

    for left, op, right, why in PROXY_RELATIONS:
        if left in document and right in document:
            lhs, rhs = document[left], document[right]
            if not (isinstance(lhs, (int, float)) and isinstance(rhs, (int, float))):
                continue
            if not math.isfinite(lhs) or not math.isfinite(rhs):
                continue
            if op == ">=" and lhs < rhs:
                problems.append(f"{left} ({lhs}) < {right} ({rhs}): {why}")

    if problems:
        raise CheckFailed("\n  ".join(["invalid proxy output:"] + problems))


def mode_proxy(args: argparse.Namespace) -> int:
    command = args.command + [
        str(args.cells),
        str(args.timesteps),
        args.kernel,
        "-f",
        "json",
    ]
    proc = run_command(command, args.timeout)
    if proc.returncode != 0:
        echo_output(proc)
        raise CheckFailed(
            f"proxy terminated unsuccessfully: {describe_exit(proc.returncode)}"
        )

    try:
        document = parse_proxy_json(extract_json_object(proc.stdout))
        check_proxy_fields(document, args.kernel)
    except CheckFailed:
        echo_output(proc)
        raise

    print(f"{args.kernel}: {len(document)} fields, all finite and in range")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="mode", required=True)

    run = sub.add_parser("run", help="run a command and assert how it terminated")
    run.add_argument(
        "--expect-failure",
        action="store_true",
        help="require unsuccessful termination (a non-zero code or a signal)",
    )
    run.add_argument(
        "--expect-output",
        default=None,
        metavar="REGEX",
        help="regular expression the combined output has to match",
    )
    run.add_argument("--timeout", type=float, default=None, help="seconds")
    run.add_argument("command", nargs=argparse.REMAINDER)
    run.set_defaults(func=mode_run)

    proxy = sub.add_parser("proxy", help="check the proxy's JSON output")
    proxy.add_argument("--kernel", required=True, help="kernel to benchmark")
    proxy.add_argument("--cells", type=int, default=10)
    proxy.add_argument("--timesteps", type=int, default=1)
    proxy.add_argument("--timeout", type=float, default=None, help="seconds")
    proxy.add_argument("command", nargs=argparse.REMAINDER)
    proxy.set_defaults(func=mode_proxy)

    args = p.parse_args(argv)
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        p.error("no command given; expected `-- <command> [args...]`")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        return args.func(args)
    except CheckFailed as error:
        print(f"smoke check failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
