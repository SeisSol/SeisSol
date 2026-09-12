#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Run one verification case and judge the result.

A case is a base parameter file plus a few overrides, run on a generated mini
mesh in its own working directory. This wrapper exists because three things
cannot be expressed in ``add_test`` alone:

* the build may not contain the configuration a case needs. The case is still
  registered so that it stays visible, and is reported as skipped rather than
  as a failure or as nothing at all.
* the parameter file has no include or override mechanism, so the overrides
  have to be applied to a copy of it here.
* a run that finishes with return code zero has not necessarily computed
  anything. The checks below insist that the output is finite and that the
  case actually did what it claims to test.

Every run leaves a small JSON record behind, which the coverage command turns
into the answer to "was every declared case executed by some build".
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# the return code CTest is told to read as "skipped"
SKIP_RETURN_CODE = 77

SECTION = re.compile(r"^\s*&(\w+)")
ASSIGNMENT = re.compile(r"^(\s*)([A-Za-z_]\w*)(\s*=\s*)(.*)$")
COMPARISON = re.compile(r"^([\w.-]+)\s*(>=|<=|!=|=)\s*(.*)$")


# --------------------------------------------------------------------------
# parameter files
# --------------------------------------------------------------------------


class ParameterFile:
    """A minimal editor for the ``&section key = value /`` parameter format.

    Comments and untouched lines are preserved verbatim, because the base file
    is meant to stay readable and reviewable. This is a stopgap: an override
    mechanism in SeisSol itself would remove the need to rewrite the file at
    all, and would also let a test state its configuration in one place
    instead of in a file plus a command line.
    """

    def __init__(self, text):
        # an ordered list of items, each either a raw line outside any section
        # or a [name, lines] section, so that rendering reproduces the input
        self.items = []
        current = None
        for line in text.splitlines():
            match = SECTION.match(line)
            if match is not None:
                current = [match.group(1), []]
                self.items.append(current)
            elif current is not None and line.strip() == "/":
                current = None
            elif current is not None:
                current[1].append(line)
            else:
                self.items.append(line)

    @property
    def sections(self):
        return [item for item in self.items if not isinstance(item, str)]

    def set(self, section, key, value):
        for name, lines in self.sections:
            if name.lower() != section.lower():
                continue
            for index, line in enumerate(lines):
                match = ASSIGNMENT.match(line)
                if match is not None and match.group(2).lower() == key.lower():
                    # keep the comment, it usually documents the units
                    comment = ""
                    position = match.group(4).find("!")
                    if position >= 0:
                        comment = " " + match.group(4)[position:]
                    lines[index] = (
                        match.group(1)
                        + match.group(2)
                        + match.group(3)
                        + value
                        + comment
                    )
                    return
            lines.append(f"{key} = {value}")
            return
        self.items.append([section, [f"{key} = {value}"]])

    def get(self, section, key):
        for name, lines in self.sections:
            if name.lower() != section.lower():
                continue
            for line in lines:
                match = ASSIGNMENT.match(line)
                if match is not None and match.group(2).lower() == key.lower():
                    value = match.group(4)
                    position = value.find("!")
                    if position >= 0:
                        value = value[:position]
                    return value.strip().strip("'\"")
        return None

    def render(self):
        parts = []
        for item in self.items:
            if isinstance(item, str):
                parts.append(item)
            else:
                parts.append(f"&{item[0]}")
                parts.extend(item[1])
                parts.append("/")
        return "\n".join(parts) + "\n"


def apply_overrides(parameters, overrides):
    for override in overrides:
        target, _, value = override.partition("=")
        if not value:
            raise SystemExit(
                f"malformed override {override!r}, expected Section.Key=Value"
            )
        section, _, key = target.rpartition(".")
        if not section:
            raise SystemExit(f"override {override!r} does not name a section")
        parameters.set(section, key, value)


# --------------------------------------------------------------------------
# capabilities
# --------------------------------------------------------------------------


def unmet_requirements(capabilities, requirements):
    """Return the requirements this build does not satisfy.

    A requirement is ``key=a|b``, ``key!=a``, ``key>=n`` or ``key<=n``. An
    unknown key is a failure of the test declaration, not of the build, so it
    raises instead of quietly skipping: a typo that skips everything would
    leave the suite green and empty.
    """
    unmet = []
    for requirement in requirements:
        match = COMPARISON.match(requirement)
        if match is None:
            raise SystemExit(f"malformed requirement {requirement!r}")
        key, operator, wanted = match.groups()
        if key not in capabilities:
            raise SystemExit(
                f"requirement {requirement!r} names an unknown capability; "
                f"known: {', '.join(sorted(capabilities))}"
            )
        have = capabilities[key]
        if operator in ("=", "!="):
            alternatives = [item.strip().lower() for item in wanted.split("|")]
            matched = str(have).lower() in alternatives
            if matched != (operator == "="):
                unmet.append(f"{key} is {have}, required {operator} {wanted}")
        else:
            try:
                left, right = float(have), float(wanted)
            except ValueError as error:
                raise SystemExit(
                    f"requirement {requirement!r} compares non-numbers"
                ) from error
            if (operator == ">=" and left < right) or (
                operator == "<=" and left > right
            ):
                unmet.append(f"{key} is {have}, required {operator} {wanted}")
    return unmet


# --------------------------------------------------------------------------
# output checks
# --------------------------------------------------------------------------


def read_receiver(path):
    """Read a receiver file into column names and rows of floats."""
    names = []
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.upper().startswith("VARIABLES"):
            names = re.findall(r'"([^"]*)"', stripped)
        elif stripped.startswith("#") or stripped.upper().startswith("TITLE"):
            continue
        else:
            rows.append([float(item) for item in stripped.split()])
    return names, rows


def check_outputs(work, prefix, require_finite, activity):
    """Check that the run produced finite output, and that it did something."""
    problems = []
    receivers = sorted(
        (work / Path(prefix).parent).glob(f"{Path(prefix).name}-receiver-*.dat")
    )
    if not receivers:
        return ["no receiver output was written"]

    extrema = {}
    for path in receivers:
        try:
            names, rows = read_receiver(path)
        except ValueError as error:
            problems.append(f"{path.name}: unreadable ({error})")
            continue
        if not rows:
            problems.append(f"{path.name}: no samples")
            continue
        for row in rows:
            for column, value in enumerate(row):
                if require_finite and not math.isfinite(value):
                    problems.append(f"{path.name}: non-finite value in column {column}")
                    break
                name = names[column] if column < len(names) else str(column)
                extrema[name] = max(extrema.get(name, 0.0), abs(value))

    for demand in activity:
        name, _, threshold = demand.partition(":")
        threshold = float(threshold)
        seen = extrema.get(name)
        if seen is None:
            problems.append(
                f"activity check wants {name!r}, which the receivers do not "
                f"contain (have: {', '.join(sorted(extrema))})"
            )
        elif seen <= threshold:
            # a case that runs cleanly but never excites anything tests nothing
            problems.append(
                f"max |{name}| is {seen:g}, expected more than {threshold:g}"
            )
    return problems


# --------------------------------------------------------------------------
# running
# --------------------------------------------------------------------------


def record(path, payload):
    if path is None:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _command_run(args):
    capabilities = json.loads(Path(args.capabilities).read_text(encoding="utf-8"))
    result = {
        "name": args.name,
        "requirements": args.require,
        "capabilities": capabilities,
    }

    unmet = unmet_requirements(capabilities, args.require)
    if unmet:
        result["status"] = "skipped"
        result["reason"] = unmet
        record(args.result, result)
        print(f"skipped: {'; '.join(unmet)}")
        return SKIP_RETURN_CODE

    work = Path(args.work_dir)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    case_dir = Path(args.case_dir)
    parameters = ParameterFile((case_dir / args.parameters).read_text(encoding="utf-8"))
    for name in args.file:
        shutil.copy(case_dir / name, work / Path(name).name)

    # an absolute mesh path keeps the working directory free of copies
    parameters.set("MeshNml", "MeshFile", f"'{Path(args.mesh).resolve()}'")
    # the contiguous distribution from the file order is the only decomposition
    # that is reproducible across rank counts and partitioner versions
    parameters.set("MeshNml", "PartitioningLib", "'none'")
    if args.timestep is not None:
        # a fixed width makes the step size independent of mesh, material and
        # order, which is what lets different builds be compared at all
        parameters.set("Discretization", "FixTimeStep", repr(args.timestep))
        if args.steps is not None:
            parameters.set("AbortCriteria", "EndTime", repr(args.steps * args.timestep))
    apply_overrides(parameters, args.set)

    prefix = parameters.get("Output", "OutputFile") or "output/out"
    (work / Path(prefix).parent).mkdir(parents=True, exist_ok=True)
    (work / "parameters.par").write_text(parameters.render(), encoding="utf-8")

    command = []
    if args.ranks > 1:
        if not args.mpiexec:
            raise SystemExit(
                "more than one rank requested but no launcher was configured"
            )
        command += [args.mpiexec, args.numproc_flag, str(args.ranks)]
        command += args.mpiexec_flags
    command += [args.binary, "parameters.par"]

    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = str(args.threads)
    environment.setdefault("OMP_PLACES", "cores")
    environment.setdefault("OMP_PROC_BIND", "close")

    print("running:", " ".join(command))
    completed = subprocess.run(command, cwd=work, env=environment, check=False)

    result["status"] = "run"
    result["command"] = command
    result["returncode"] = completed.returncode
    problems = []
    if completed.returncode != 0:
        problems.append(f"exit code {completed.returncode}")
    else:
        problems += check_outputs(work, prefix, args.require_finite, args.activity)

    result["problems"] = problems
    record(args.result, result)
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    return 1 if problems else 0


def _command_coverage(args):
    """Check that every declared case was run somewhere, not just skipped."""
    declared = json.loads(Path(args.declared).read_text(encoding="utf-8"))
    results = {}
    for path in sorted(Path(args.results).glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        results[payload["name"]] = payload

    ran = [name for name, payload in results.items() if payload["status"] == "run"]
    skipped = {
        name: payload["reason"]
        for name, payload in results.items()
        if payload["status"] == "skipped"
    }
    absent = [case["name"] for case in declared["cases"] if case["name"] not in results]

    print(f"declared {len(declared['cases'])}, ran {len(ran)}, skipped {len(skipped)}")
    for name, reason in sorted(skipped.items()):
        print(f"  skipped {name}: {'; '.join(reason)}")
    for name in absent:
        print(f"  no result {name}")

    if args.history:
        history = Path(args.history)
        seen = (
            json.loads(history.read_text(encoding="utf-8")) if history.exists() else {}
        )
        for name in ran:
            seen[name] = True
        history.write_text(json.dumps(seen, indent=2, sort_keys=True) + "\n", "utf-8")
        never = [
            case["name"] for case in declared["cases"] if not seen.get(case["name"])
        ]
        if never:
            # a case no build has ever run is a hole in the matrix, and it looks
            # exactly like a passing suite unless it is called out
            print("never run by any build:", file=sys.stderr)
            for name in sorted(never):
                print(f"  {name}", file=sys.stderr)
            return 1

    if absent and args.require_results:
        return 1
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="run a SeisSol verification case")
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="run a single case")
    run.add_argument("--name", required=True)
    run.add_argument("--binary", required=True)
    run.add_argument("--capabilities", required=True)
    run.add_argument("--case-dir", required=True)
    run.add_argument("--work-dir", required=True)
    run.add_argument("--mesh", required=True)
    run.add_argument("--parameters", required=True)
    run.add_argument("--result")
    run.add_argument("--ranks", type=int, default=1)
    run.add_argument("--threads", type=int, default=1)
    run.add_argument("--steps", type=int)
    run.add_argument("--timestep", type=float)
    run.add_argument("--mpiexec", default="")
    run.add_argument("--numproc-flag", default="-n")
    run.add_argument("--mpiexec-flags", nargs="*", default=[])
    run.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE")
    run.add_argument("--file", action="append", default=[], metavar="NAME")
    run.add_argument("--require", action="append", default=[], metavar="KEY=VALUE")
    run.add_argument(
        "--activity",
        action="append",
        default=[],
        metavar="COLUMN:THRESHOLD",
        help="insist that the run actually excited something",
    )
    run.add_argument("--no-finite-check", dest="require_finite", action="store_false")
    run.set_defaults(require_finite=True, func=_command_run)

    coverage = commands.add_parser(
        "coverage", help="account for declared versus run cases"
    )
    coverage.add_argument("--declared", required=True)
    coverage.add_argument("--results", required=True)
    coverage.add_argument(
        "--history", help="accumulate across builds and fail on holes"
    )
    coverage.add_argument("--require-results", action="store_true")
    coverage.set_defaults(func=_command_coverage)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
