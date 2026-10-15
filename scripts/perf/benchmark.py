#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""Measure the kernels of SeisSol with its proxy, and compare two measurements.

``run``
    Runs the proxy (``proxyseissol``) for each configuration and kernel asked
    for, a few times after a warm-up, and writes every sample to a JSON
    report, together with what the machine looked like while it ran (CPU,
    threads, frequency governor, load). The time per cell and step is the
    measure of a benchmark; the GFLOP/s and GiB/s the proxy reports go along.

``compare``
    Holds a report against a baseline, benchmark by benchmark, on the median
    time. A change counts as a regression if it exceeds both the threshold
    and the noise the two reports measured themselves (the spread of their
    repetitions), so that a noisy machine does not report a regression on
    every run. Reports from machines that look different are not compared.

``noise``
    Prints the spread of the repetitions of a report, which tells whether a
    machine measures steadily enough to compare runs on it at all, and which
    threshold it allows.

The proxy is a deterministic workload without I/O and communication, which is
what makes runs comparable over time. Everything else that a simulation spends
its time in is not covered here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import shlex
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPORT_SCHEMA = "seissol-benchmark/1"

# the kernels of the proxy that together make up a time step of a simulation without and with
# dynamic rupture; the others are parts of these
DEFAULT_KERNELS = ("all", "all_dr")

# what the reports of two machines have to agree in to be compared
COMPARABLE_KEYS = ("cpu", "cpus", "omp_num_threads", "omp_places", "omp_proc_bind")

# a configuration name as describeBuiltConfigs() starts it, at the start of an unindented line
CONFIG_NAME = re.compile(
    r"^([a-z]+-[a-z]+(?:-m\d+)?-o\d+-f(?:32|64)-[a-z]+(?:-s\d+)?)\s*$"
)


class BenchmarkError(Exception):
    """A run or a comparison could not be done."""


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def median(values: list[float]) -> float:
    return statistics.median(values)


def robust_spread(values: list[float]) -> float:
    """The relative spread of the values: 1.4826 * MAD / median (the CV for normal noise).

    The median absolute deviation ignores the odd outlier, as one repetition that another process
    on the machine slowed down; with a single value, there is no spread to measure.
    """
    if len(values) < 2:
        return 0.0
    center = median(values)
    if center == 0:
        return 0.0
    deviation = median([abs(value - center) for value in values])
    return 1.4826 * deviation / center


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------


def _read(path: str) -> str | None:
    try:
        return Path(path).read_text(encoding="utf-8").strip()
    except OSError:
        return None


def cpu_model() -> str:
    text = _read("/proc/cpuinfo") or ""
    for key in ("model name", "Model", "cpu model", "Hardware"):
        match = re.search(rf"^{key}\s*:\s*(.+)$", text, re.MULTILINE)
        if match:
            return match.group(1).strip()
    return platform.processor() or platform.machine()


def describe_environment() -> dict[str, Any]:
    """What a measurement depends on besides the code: the machine, and how it ran."""
    boost = _read("/sys/devices/system/cpu/cpufreq/boost")
    no_turbo = _read("/sys/devices/system/cpu/intel_pstate/no_turbo")
    if boost is not None:
        turbo: str | None = "on" if boost == "1" else "off"
    elif no_turbo is not None:
        turbo = "off" if no_turbo == "1" else "on"
    else:
        turbo = None
    return {
        "cpu": cpu_model(),
        "cpus": os.cpu_count(),
        "machine": platform.machine(),
        "kernel": platform.release(),
        "governor": _read("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "turbo": turbo,
        "load": list(os.getloadavg()) if hasattr(os, "getloadavg") else None,
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "omp_places": os.environ.get("OMP_PLACES"),
        "omp_proc_bind": os.environ.get("OMP_PROC_BIND"),
        "commit": os.environ.get("GITHUB_SHA"),
    }


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# Running the proxy
# ---------------------------------------------------------------------------


def built_configs(proxy: list[str]) -> list[str]:
    """The configurations the proxy holds, as its help lists them."""
    completed = subprocess.run(
        proxy + ["--help"], capture_output=True, text=True, check=False
    )
    names = []
    for line in (completed.stdout + completed.stderr).splitlines():
        match = CONFIG_NAME.match(line)
        if match and match.group(1) not in names:
            names.append(match.group(1))
    if not names:
        raise BenchmarkError(
            f"{shlex.join(proxy)} --help lists no configuration; name them with --config"
        )
    return names


def parse_proxy_output(stdout: str) -> dict[str, Any]:
    """The JSON object the proxy writes as the last line of its output in JSON mode."""
    lines = [line for line in stdout.splitlines() if line.strip()]
    if not lines or not lines[-1].lstrip().startswith("{"):
        raise BenchmarkError("the proxy wrote no JSON object")
    return json.loads(lines[-1])


def run_proxy(
    proxy: list[str], config: str, kernel: str, cells: int, timesteps: int
) -> dict[str, Any]:
    command = proxy + [str(cells), str(timesteps), kernel, "-f", "json", "-c", config]
    start = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        raise BenchmarkError(
            f"{shlex.join(command)} ended with exit code {completed.returncode}:\n"
            + completed.stdout[-2000:]
            + completed.stderr[-2000:]
        )
    document = parse_proxy_output(completed.stdout)
    document["wall"] = elapsed
    return document


def benchmark_name(config: str, kernel: str) -> str:
    return f"{config}/{kernel}"


def mode_run(args: argparse.Namespace) -> int:
    # absolute, as a bare name (also ./proxyseissol as a Path) would be looked up in the PATH
    proxy = shlex.split(args.launcher) + [str(args.proxy.resolve())]
    configs = args.config or built_configs(proxy)
    kernels = args.kernel or list(DEFAULT_KERNELS)
    report: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        "environment": describe_environment(),
        "proxy": {"path": str(args.proxy), "sha256": file_digest(args.proxy)},
        "parameters": {
            "cells": args.cells,
            "timesteps": args.timesteps,
            "repetitions": args.repetitions,
            "warmup": args.warmup,
        },
        "benchmarks": {},
    }
    cell_steps = args.cells * args.timesteps
    for config in configs:
        for kernel in kernels:
            for _ in range(args.warmup):
                run_proxy(proxy, config, kernel, args.cells, args.timesteps)
            samples = [
                run_proxy(proxy, config, kernel, args.cells, args.timesteps)
                for _ in range(args.repetitions)
            ]
            times = [sample["time"] / cell_steps * 1e9 for sample in samples]
            name = benchmark_name(config, kernel)
            report["benchmarks"][name] = {
                "config": config,
                "kernel": kernel,
                "ns_per_cell_step": times,
                "gflops_hw": [sample["gflops-hw"] for sample in samples],
                "gflops_nz": [sample["gflops-nz"] for sample in samples],
                "gibs": [sample["gibs"] for sample in samples],
            }
            print(
                f"{name}: {median(times):.1f} ns/(cell step), "
                f"spread {100 * robust_spread(times):.2f} %, "
                f"{median(report['benchmarks'][name]['gflops_hw']):.1f} HW-GFLOP/s",
                flush=True,
            )
    report["environment"]["load_after"] = (
        list(os.getloadavg()) if hasattr(os, "getloadavg") else None
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.github_benchmark:
        args.github_benchmark.write_text(
            json.dumps(github_benchmark_entries(report), indent=2) + "\n",
            encoding="utf-8",
        )
    return 0


def github_benchmark_entries(report: dict[str, Any]) -> list[dict[str, Any]]:
    """The report in the customSmallerIsBetter format of github-action-benchmark."""
    entries = []
    for name, data in report["benchmarks"].items():
        times = data["ns_per_cell_step"]
        entries.append(
            {
                "name": name,
                "unit": "ns/(cell step)",
                "value": round(median(times), 3),
                "range": f"{100 * robust_spread(times):.2f}%",
                "extra": f"{median(data['gflops_hw']):.2f} HW-GFLOP/s",
            }
        )
    return entries


# ---------------------------------------------------------------------------
# Comparing two reports
# ---------------------------------------------------------------------------


def load_report(path: Path) -> dict[str, Any]:
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("schema") != REPORT_SCHEMA:
        raise BenchmarkError(f"{path} is no report of the schema {REPORT_SCHEMA}")
    return report


def environment_differences(
    baseline: dict[str, Any], current: dict[str, Any]
) -> list[str]:
    return [
        f"{key}: {baseline['environment'].get(key)} -> {current['environment'].get(key)}"
        for key in COMPARABLE_KEYS
        if baseline["environment"].get(key) != current["environment"].get(key)
    ]


def compare_reports(
    baseline: dict[str, Any],
    current: dict[str, Any],
    threshold: float,
    noise_factor: float,
) -> list[dict[str, Any]]:
    """Per benchmark of both reports: the change of the median time and its verdict.

    A change is significant if it exceeds the threshold as well as noise_factor times the
    combined spread of the two reports' repetitions.
    """
    rows = []
    for name in sorted(set(baseline["benchmarks"]) | set(current["benchmarks"])):
        before = baseline["benchmarks"].get(name)
        after = current["benchmarks"].get(name)
        if before is None or after is None:
            rows.append(
                {
                    "name": name,
                    "verdict": "new" if before is None else "gone",
                }
            )
            continue
        old = median(before["ns_per_cell_step"])
        new = median(after["ns_per_cell_step"])
        change = new / old - 1.0
        noise = math.hypot(
            robust_spread(before["ns_per_cell_step"]),
            robust_spread(after["ns_per_cell_step"]),
        )
        limit = max(threshold, noise_factor * noise)
        if change > limit:
            verdict = "slower"
        elif change < -limit:
            verdict = "faster"
        else:
            verdict = "unchanged"
        rows.append(
            {
                "name": name,
                "before": old,
                "after": new,
                "change": change,
                "limit": limit,
                "verdict": verdict,
            }
        )
    return rows


def render_comparison(rows: list[dict[str, Any]], differences: list[str]) -> str:
    lines = ["### Kernel performance (proxy)", ""]
    if differences:
        lines += [
            "The reports come from different machines or settings and were not compared:",
            "",
        ]
        lines += [f"- {difference}" for difference in differences]
        return "\n".join(lines) + "\n"
    lines += [
        "| benchmark | baseline [ns/(cell step)] | now | change | limit | |",
        "|---|---:|---:|---:|---:|---|",
    ]
    marks = {"slower": "🔴 slower", "faster": "🟢 faster", "unchanged": ""}
    for row in rows:
        if "change" not in row:
            lines.append(f"| {row['name']} | | | | | {row['verdict']} |")
            continue
        lines.append(
            f"| {row['name']} | {row['before']:.1f} | {row['after']:.1f} "
            f"| {100 * row['change']:+.2f} % | ±{100 * row['limit']:.2f} % "
            f"| {marks[row['verdict']]} |"
        )
    return "\n".join(lines) + "\n"


def mode_compare(args: argparse.Namespace) -> int:
    baseline = load_report(args.baseline)
    current = load_report(args.current)
    differences = environment_differences(baseline, current)
    rows = (
        []
        if differences
        else compare_reports(baseline, current, args.threshold, args.noise_factor)
    )
    text = render_comparison(rows, differences)
    sys.stdout.write(text)
    if args.summary:
        with args.summary.open("a", encoding="utf-8") as handle:
            handle.write(text)
    slower = [row["name"] for row in rows if row["verdict"] == "slower"]
    if slower and args.fail_on_regression:
        print(f"slower: {', '.join(slower)}", file=sys.stderr)
        return 1
    return 0


def mode_noise(args: argparse.Namespace) -> int:
    report = load_report(args.report)
    print(f"{'benchmark':<55} {'median':>10} {'spread':>8} {'min..max':>18}")
    worst = 0.0
    for name, data in sorted(report["benchmarks"].items()):
        times = data["ns_per_cell_step"]
        spread = robust_spread(times)
        worst = max(worst, spread)
        span = f"{min(times):.1f}..{max(times):.1f}"
        print(f"{name:<55} {median(times):>10.1f} {100 * spread:>7.2f}% {span:>18}")
    print(
        f"largest spread: {100 * worst:.2f} %; a comparison with a noise factor of "
        f"{args.noise_factor:g} tells changes above "
        f"{100 * args.noise_factor * math.sqrt(2) * worst:.2f} % apart"
    )
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="mode", required=True)

    run = sub.add_parser("run", help="measure the kernels with the proxy")
    run.add_argument("--proxy", type=Path, required=True, help="the proxy executable")
    run.add_argument(
        "--launcher",
        default="",
        help="command to start the proxy with, e.g. 'numactl --cpunodebind=0'",
    )
    run.add_argument(
        "--config",
        action="append",
        help="configuration to measure (repeatable; default: all of the proxy)",
    )
    run.add_argument(
        "--kernel",
        action="append",
        help=f"kernel to measure (repeatable; default: {', '.join(DEFAULT_KERNELS)})",
    )
    run.add_argument("--cells", type=int, default=50000)
    run.add_argument("--timesteps", type=int, default=20)
    run.add_argument("--repetitions", type=int, default=5)
    run.add_argument("--warmup", type=int, default=1)
    run.add_argument("--output", type=Path, required=True, help="report to write")
    run.add_argument(
        "--github-benchmark",
        type=Path,
        help="also write the medians for github-action-benchmark (customSmallerIsBetter)",
    )
    run.set_defaults(func=mode_run)

    compare = sub.add_parser("compare", help="compare a report with a baseline")
    compare.add_argument("--baseline", type=Path, required=True)
    compare.add_argument("--current", type=Path, required=True)
    compare.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="smallest relative change that counts (default: 0.05)",
    )
    compare.add_argument(
        "--noise-factor",
        type=float,
        default=3.0,
        help="how many times the combined spread a change has to exceed (default: 3)",
    )
    compare.add_argument("--summary", type=Path, help="append the table to this file")
    compare.add_argument("--fail-on-regression", action="store_true")
    compare.set_defaults(func=mode_compare)

    noise = sub.add_parser("noise", help="print the spread of the repetitions")
    noise.add_argument("report", type=Path)
    noise.add_argument("--noise-factor", type=float, default=3.0)
    noise.set_defaults(func=mode_noise)

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        return args.func(args)
    except BenchmarkError as error:
        print(f"benchmark: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
