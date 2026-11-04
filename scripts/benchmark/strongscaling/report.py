# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The scaling of a campaign: tables, plots and a report, from results/runs.csv.

A series is a build and a variant over the node counts. Per node count, its
time is the median over the repetitions that completed, with their range; its
speedup and parallel efficiency are taken against the smallest node count of
the series, and its throughput is the cell updates per second and node.

    results/scaling.csv   a row per series and node count
    results/scaling.json  the same, with the checks of the results
    results/report.md     the tables, the checks and the plots
    results/plots/        time, efficiency, throughput, breakdown, imbalance
"""

from __future__ import annotations

import csv
import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .campaign import COMPLETED
from .campaign import Layout as Directory
from .campaign import utc_now
from .collect import number, write_csv
from .config import Campaign
from .plan import Plan

METRICS = {
    "total": "the simulation time, of its slowest rank",
    "compute": "the compute part of the simulation time, of its slowest rank",
    "steady": "the epochs after the warm-up, scaled to the whole simulation",
}

SCALING_COLUMNS = [
    "build",
    "variant",
    "nodes",
    "ranks",
    "cellsPerRank",
    "runs",
    "completed",
    "time",
    "timeMin",
    "timeMax",
    "speedup",
    "efficiency",
    "cellUpdatesPerSecondPerNode",
    "duplicateShare",
    "loadImbalance",
    "outsideKernelsShare",
    "ioShare",
    "initTime",
    "hostMemoryPeak",
]


def load_runs(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key, value in list(row.items()):
            if key in (
                "run",
                "build",
                "variant",
                "status",
                "runId",
                "job",
                "version",
                "commit",
                "problems",
            ):
                continue
            converted = number(value)
            row[key] = converted if converted is not None else (value or None)
    return rows


def metric_value(row: dict[str, Any], metric: str) -> float | None:
    if metric == "total":
        return row.get("wallTime")
    if metric == "compute":
        return row.get("computeTime")
    wall = row.get("steadyWallTime")
    simulated = row.get("steadySimulatedTime")
    end = row.get("endTime")
    if not wall or not simulated or not end:
        return None
    return wall / simulated * end


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


@dataclass
class Point:
    nodes: int
    ranks: int
    runs: int
    times: list[float]
    rows: list[dict[str, Any]]
    time: float | None = None
    speedup: float | None = None
    efficiency: float | None = None

    def median_of(self, key: str) -> float | None:
        return median(
            [row[key] for row in self.rows if isinstance(row.get(key), float)]
        )


@dataclass
class Series:
    build: str
    variant: str
    points: list[Point] = field(default_factory=list)

    @property
    def label(self) -> str:
        return f"{self.build}/{self.variant}"


def make_series(rows: list[dict[str, Any]], metric: str) -> list[Series]:
    order: list[tuple[str, str]] = []
    for row in rows:
        key = (row["build"], row["variant"])
        if key not in order:
            order.append(key)
    result = []
    for build, variant in order:
        series = Series(build, variant)
        mine = [r for r in rows if r["build"] == build and r["variant"] == variant]
        for nodes in sorted({int(r["nodes"]) for r in mine}):
            at = [r for r in mine if int(r["nodes"]) == nodes]
            done = [r for r in at if r["status"] == COMPLETED]
            times = [t for t in (metric_value(r, metric) for r in done) if t]
            series.points.append(
                Point(
                    nodes=nodes,
                    ranks=int(at[0]["ranks"]),
                    runs=len(at),
                    times=times,
                    rows=done,
                    time=median(times),
                )
            )
        base = next((p for p in series.points if p.time), None)
        if base is not None:
            for point in series.points:
                if point.time:
                    point.speedup = base.time / point.time
                    point.efficiency = point.speedup * base.nodes / point.nodes
        result.append(series)
    return result


def scaling_rows(all_series: list[Series]) -> list[dict[str, Any]]:
    rows = []
    for series in all_series:
        for point in series.points:
            updates = point.median_of("cellUpdates")
            duplicates = point.median_of("duplicateCellUpdates")
            outside = point.median_of("outsideKernelsMax")
            compute = point.median_of("computeTime")
            io = point.median_of("ioTime")
            wall = point.median_of("wallTime")
            rows.append(
                {
                    "build": series.build,
                    "variant": series.variant,
                    "nodes": point.nodes,
                    "ranks": point.ranks,
                    "cellsPerRank": point.median_of("cellsPerRank"),
                    "runs": point.runs,
                    "completed": len(point.rows),
                    "time": point.time,
                    "timeMin": min(point.times) if point.times else None,
                    "timeMax": max(point.times) if point.times else None,
                    "speedup": point.speedup,
                    "efficiency": point.efficiency,
                    "cellUpdatesPerSecondPerNode": (
                        updates / point.time / point.nodes
                        if updates and point.time
                        else None
                    ),
                    "duplicateShare": (
                        duplicates / updates
                        if duplicates is not None and updates
                        else None
                    ),
                    "loadImbalance": point.median_of("loadImbalance"),
                    "outsideKernelsShare": (
                        outside / compute if outside and compute else None
                    ),
                    "ioShare": io / wall if io is not None and wall else None,
                    "initTime": point.median_of("initTime"),
                    "hostMemoryPeak": point.median_of("hostMemoryPeak"),
                }
            )
    return rows


def checks(rows: list[dict[str, Any]], all_series: list[Series]) -> list[str]:
    """What looks wrong in the results."""
    findings = []
    for row in rows:
        if row["status"] != COMPLETED:
            findings.append(
                f"{row['run']}: {row['status']}"
                + (f" ({row['problems']})" if row.get("problems") else "")
            )
        elif row.get("problems"):
            findings.append(f"{row['run']}: {row['problems']}")
    for series in all_series:
        updates = {
            point.nodes: point.median_of("cellUpdates")
            for point in series.points
            if point.rows
        }
        values = {v for v in updates.values() if v}
        if len(values) > 1:
            findings.append(
                f"{series.label}: the cell updates differ between the node counts "
                f"({', '.join(f'{n}: {u:.0f}' for n, u in sorted(updates.items()) if u)}); "
                "the runs did not compute the same, e.g. not to the same end time"
            )
        for point in series.points:
            if len(point.times) > 1 and point.time:
                spread = (max(point.times) - min(point.times)) / point.time
                if spread > 0.1:
                    findings.append(
                        f"{series.label} on {point.nodes} nodes: the repetitions spread by "
                        f"{100 * spread:.0f}% of their median"
                    )
    return findings


def fmt(value: Any, digits: int = 3) -> str:
    if value is None or value == "":
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def plots(
    all_series: list[Series], table: list[dict[str, Any]], folder: Path, metric: str
) -> list[Path]:
    """The plots of the scaling; none if matplotlib is not there."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    drawn = [s for s in all_series if any(p.time for p in s.points)]
    if not drawn:
        return []

    node_counts = sorted({p.nodes for s in drawn for p in s.points})

    def axes(title: str, ylabel: str, logy: bool = False):
        figure, axis = plt.subplots(figsize=(6.4, 4.4))
        axis.set_xscale("log", base=2)
        if logy:
            axis.set_yscale("log")
        # the node counts as they are, rather than as powers of two
        axis.set_xticks(node_counts)
        axis.set_xticklabels([str(n) for n in node_counts])
        axis.minorticks_off()
        axis.set_xlabel("nodes")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(True, which="both", alpha=0.3)
        return figure, axis

    def save(figure, name: str) -> None:
        path = folder / name
        figure.tight_layout()
        figure.savefig(path, dpi=150)
        plt.close(figure)
        written.append(path)

    figure, axis = axes(f"Time ({metric})", "time [s]", logy=True)
    for series in drawn:
        points = [p for p in series.points if p.time]
        nodes = [p.nodes for p in points]
        times = [p.time for p in points]
        low = [p.time - min(p.times) for p in points]
        high = [max(p.times) - p.time for p in points]
        line = axis.errorbar(
            nodes, times, yerr=[low, high], marker="o", capsize=3, label=series.label
        )
        ideal = [points[0].time * points[0].nodes / n for n in nodes]
        axis.plot(nodes, ideal, linestyle="--", color=line[0].get_color(), alpha=0.5)
    axis.legend()
    save(figure, "time.png")

    figure, axis = axes("Parallel efficiency", "efficiency")
    for series in drawn:
        points = [p for p in series.points if p.efficiency]
        axis.plot(
            [p.nodes for p in points],
            [p.efficiency for p in points],
            marker="o",
            label=series.label,
        )
    axis.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
    axis.set_ylim(bottom=0)
    axis.legend()
    save(figure, "efficiency.png")

    throughput = [r for r in table if r.get("cellUpdatesPerSecondPerNode")]
    if throughput:
        figure, axis = axes("Throughput", "cell updates per second and node")
        for series in drawn:
            mine = [
                r
                for r in throughput
                if r["build"] == series.build and r["variant"] == series.variant
            ]
            axis.plot(
                [r["nodes"] for r in mine],
                [r["cellUpdatesPerSecondPerNode"] for r in mine],
                marker="o",
                label=series.label,
            )
        axis.set_ylim(bottom=0)
        axis.legend()
        save(figure, "throughput.png")

    imbalance = [r for r in table if r.get("loadImbalance") is not None]
    if imbalance:
        figure, axis = axes("Load imbalance of the kernels", "1 - mean / max")
        for series in drawn:
            mine = [
                r
                for r in imbalance
                if r["build"] == series.build and r["variant"] == series.variant
            ]
            axis.plot(
                [r["nodes"] for r in mine],
                [r["loadImbalance"] for r in mine],
                marker="o",
                label=series.label,
            )
        axis.set_ylim(bottom=0)
        axis.legend()
        save(figure, "imbalance.png")

    # where the time of a rank goes: the kernels, the rest of the compute time, the blocking IO
    count = len(drawn)
    figure, all_axes = plt.subplots(
        1, count, figsize=(4.0 * count + 1.0, 4.2), squeeze=False
    )
    for axis, series in zip(all_axes[0], drawn):
        points = [p for p in series.points if p.rows]
        labels = [str(p.nodes) for p in points]
        parts = [
            [
                p.median_of("kernelTimeMean") or 0.0,
                p.median_of("outsideKernelsMean") or 0.0,
                p.median_of("ioTime") or 0.0,
            ]
            for p in points
        ]
        shares = [[100 * v / (sum(part) or 1.0) for v in part] for part in parts]
        bottom = [0.0] * len(points)
        for index, label in enumerate(
            ("kernels", "outside of the kernels", "blocking IO")
        ):
            values = [share[index] for share in shares]
            axis.bar(labels, values, bottom=bottom, label=label)
            bottom = [b + v for b, v in zip(bottom, values)]
        axis.set_ylim(0, 100)
        axis.set_xlabel("nodes")
        axis.set_ylabel("share of the time of a rank [%]")
        axis.set_title(series.label)
    all_axes[0][0].legend(loc="lower left")
    save(figure, "breakdown.png")
    return written


def report(
    campaign: Campaign, plan: Plan, directory: Directory, metric: str = "total"
) -> list[Path]:
    """Writes the scaling tables, plots and report of the campaign; returns what it wrote."""
    results = directory.results
    rows = load_runs(results / "runs.csv")
    all_series = make_series(rows, metric)
    table = scaling_rows(all_series)
    findings = checks(rows, all_series)
    written = [write_csv(results / "scaling.csv", SCALING_COLUMNS, table)]
    (results / "scaling.json").write_text(
        json.dumps(
            {
                "schema": "seissol-scaling-report/1",
                "campaign": campaign.name,
                "created": utc_now(),
                "metric": metric,
                "series": table,
                "checks": findings,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    written.append(results / "scaling.json")
    plot_paths = plots(all_series, table, results / "plots", metric)
    written += plot_paths

    commits = sorted({str(r["commit"]) for r in rows if r.get("commit")})
    versions = sorted({str(r["version"]) for r in rows if r.get("version")})
    lines = [
        f"# Strong scaling: {campaign.name}",
        "",
        f"- machine: {campaign.machine.name} ({campaign.machine.scheduler})",
        f"- case: {campaign.case.parameter_file}",
        f"- cells: {fmt(plan.cells)}",
        f"- SeisSol: {', '.join(versions) or '-'} ({', '.join(c[:12] for c in commits) or '-'})",
        f"- time: {METRICS[metric]}; per node count the median over the repetitions, "
        "with their range",
        "- cell updates: one per cell and timestep; duplicates: the further updates of copy "
        "cells, once per further neighbor, relative to them; outside kernels: the compute "
        "time outside of the compute kernels (waiting, communication), relative to it; IO: "
        "the blocking IO, relative to the simulation time",
        f"- written: {utc_now()}",
        "",
    ]
    for series in all_series:
        lines += [
            f"## {series.label}",
            "",
            "| nodes | ranks | cells/rank | runs | time [s] | range [s] | speedup | efficiency | "
            "cell updates/s/node | duplicates | imbalance | outside kernels | IO |",
            "|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for entry in table:
            if entry["build"] != series.build or entry["variant"] != series.variant:
                continue
            spread = (
                f"{fmt(entry['timeMin'])} – {fmt(entry['timeMax'])}"
                if entry["timeMin"] is not None
                else "-"
            )
            share = entry["outsideKernelsShare"]
            io = entry["ioShare"]
            duplicates = entry["duplicateShare"]
            lines.append(
                f"| {entry['nodes']} | {entry['ranks']} | {fmt(entry['cellsPerRank'])} | "
                f"{entry['completed']}/{entry['runs']} | {fmt(entry['time'])} | {spread} | "
                f"{fmt(entry['speedup'])} | {fmt(entry['efficiency'])} | "
                f"{fmt(entry['cellUpdatesPerSecondPerNode'])} | "
                f"{fmt(100 * duplicates if duplicates is not None else None, 2)}% | "
                f"{fmt(entry['loadImbalance'])} | "
                f"{fmt(100 * share if share is not None else None, 2)}% | "
                f"{fmt(100 * io if io is not None else None, 2)}% |"
            )
        lines.append("")
    lines += ["## Checks", ""]
    lines += [f"- {finding}" for finding in findings] or ["- nothing stands out"]
    lines.append("")
    if plot_paths:
        lines += ["## Plots", ""]
        lines += [f"![{p.stem}](plots/{p.name})" for p in plot_paths]
    else:
        lines += ["(no plots: matplotlib is not installed)"]
    lines.append("")
    path = results / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    written.append(path)
    return written
