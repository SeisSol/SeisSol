# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""What the runs of a campaign wrote, gathered into tables.

    results/runs.csv     a row per run: its setting, how it went, its times
    results/ranks.csv    a row per rank of every run that completed
    results/epochs.csv   a row per epoch of every run that completed

A run is taken from its latest attempt that completed, or from its latest
attempt if none did. Its run ID ties the files of an attempt together: a file
of another run ID is not taken.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from .campaign import COMPLETED
from .campaign import Layout as Directory
from .campaign import State, attempt_status, read_json
from .config import Campaign
from .plan import Plan, Run

RUN_COLUMNS = [
    "run",
    "build",
    "variant",
    "repetition",
    "nodes",
    "ranks",
    "ranksPerNode",
    "threadsPerRank",
    "attempt",
    "attempts",
    "runId",
    "job",
    "status",
    "hosts",
    "cells",
    "cellsPerRank",
    "endTime",
    "currentTime",
    "initTime",
    "wallTime",
    "computeTime",
    "ioTime",
    "kernelTimeMax",
    "kernelTimeMean",
    "outsideKernelsMax",
    "outsideKernelsMean",
    "loadImbalance",
    "steadyEpochs",
    "steadyWallTime",
    "steadySimulatedTime",
    "cellUpdates",
    "duplicateCellUpdates",
    "cellUpdatesPerSecond",
    "hardwareFlopPerSecond",
    "nonzeroFlopPerSecond",
    "hostMemoryPeak",
    "version",
    "commit",
    "pinningProblems",
    "problems",
]


def get(data: Any, path: str, default: Any = None) -> Any:
    """The member at the dotted @p path of nested mappings, or @p default."""
    value = data
    for key in path.split("."):
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))
    except OSError:
        return []


def write_csv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: "" if row.get(key) is None else row[key] for key in columns}
            )
    return path


def number(text: Any) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def cpu_set(text: str) -> set[int] | None:
    """The CPUs of a mask as the pinning table writes it, e.g. 0,1,2; None if unknown."""
    text = (text or "").strip()
    if not text or text in ("none", "no-info"):
        return None
    cpus = set()
    for part in text.split(","):
        part = part.strip()
        if "-" in part:
            low, high = part.split("-", 1)
            cpus.update(range(int(low), int(high) + 1))
        elif part:
            cpus.add(int(part))
    return cpus


def pinning_problems(rows: list[dict[str, str]], threads_per_rank: int) -> list[str]:
    """What is wrong with the pinning of a run, from its threadPinning table."""
    problems = []
    by_host: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_host.setdefault(row.get("hostname", ""), []).append(row)
    for host, ranks in sorted(by_host.items()):
        workers = {}
        for row in ranks:
            rank = row.get("rank", "?")
            cpus = cpu_set(row.get("workermask", ""))
            if cpus is None:
                continue
            workers[rank] = cpus
            if len(cpus) < threads_per_rank:
                problems.append(
                    f"{host}: rank {rank} has {len(cpus)} CPUs for {threads_per_rank} threads"
                )
            numa = {
                n
                for n in (row.get("workernuma") or "").split(",")
                if n.strip().isdigit()
            }
            if len(numa) > 1:
                problems.append(
                    f"{host}: the workers of rank {rank} span NUMA domains {','.join(sorted(numa))}"
                )
        names = sorted(workers, key=lambda r: int(r) if r.isdigit() else 0)
        for i, first in enumerate(names):
            for second in names[i + 1 :]:
                shared = workers[first] & workers[second]
                if shared:
                    problems.append(
                        f"{host}: ranks {first} and {second} share {len(shared)} CPU(s)"
                    )
        for row in ranks:
            comm = cpu_set(row.get("commthread_mask", ""))
            if not comm:
                continue
            for rank, cpus in workers.items():
                if comm & cpus:
                    problems.append(
                        f"{host}: the communication thread of rank {row.get('rank')} shares "
                        f"CPUs with the workers of rank {rank}"
                    )
    return problems


def steady_state(
    epochs: list[dict[str, str]], warmup: int
) -> tuple[int, float | None, float | None]:
    """The epochs after the first @p warmup ones: their count, wall and simulated time."""
    rows = sorted(epochs, key=lambda row: number(row.get("epoch")) or 0)[warmup:]
    if not rows:
        return 0, None, None
    wall = sum(number(row.get("wallMax")) or 0.0 for row in rows)
    simulated = sum(
        (number(row.get("timeTo")) or 0.0) - (number(row.get("timeFrom")) or 0.0)
        for row in rows
    )
    return len(rows), wall, simulated


def chosen_attempt(
    run: Run, directory: Directory, state: State, jobs: dict[str, str]
) -> tuple[dict[str, Any] | None, str, int]:
    attempts = [
        a for a in state.attempts(run.name) if a["definition"] == run.definition
    ]
    latest = None
    for attempt in reversed(attempts):
        status = attempt_status(
            directory.attempt(run.name, attempt["attempt"]),
            jobs.get(attempt.get("job")),
        )
        if latest is None:
            latest = (attempt, status)
        if status == COMPLETED:
            return attempt, status, len(attempts)
    if latest is None:
        return None, "planned", 0
    return latest[0], latest[1], len(attempts)


def run_row(
    campaign: Campaign,
    run: Run,
    attempt: dict[str, Any] | None,
    status: str,
    count: int,
    directory: Directory,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    row: dict[str, Any] = {
        "run": run.name,
        "build": run.build,
        "variant": run.variant,
        "repetition": run.repetition,
        "nodes": run.nodes,
        "ranks": run.ranks,
        "ranksPerNode": run.layout.ranks_per_node,
        "threadsPerRank": run.layout.threads_per_rank,
        "attempt": attempt["attempt"] if attempt else None,
        "attempts": count,
        "runId": attempt["runId"] if attempt else None,
        "job": attempt.get("job") if attempt else None,
        "status": status,
    }
    if attempt is None:
        return row, [], []
    folder = directory.attempt(run.name, attempt["attempt"])
    output = folder / "output"
    problems = []
    try:
        summary = read_json(output / "data-summary.json")
    except (OSError, ValueError):
        summary = None
    if summary is None:
        return row, [], []
    if summary.get("runId") != attempt["runId"]:
        row["problems"] = (
            f"the summary is of run ID {summary.get('runId')}, not of this attempt"
        )
        return row, [], []
    row.update(
        {
            "hosts": get(summary, "parallel.nodes"),
            "cells": get(summary, "mesh.cells.sum"),
            "endTime": get(summary, "simulation.endTime"),
            "currentTime": get(summary, "simulation.currentTime"),
            "initTime": get(summary, "time.init.total.max"),
            "wallTime": get(summary, "time.simulation.total.max"),
            "computeTime": get(summary, "time.simulation.compute.max"),
            "ioTime": get(summary, "time.simulation.io.max"),
            "kernelTimeMax": get(summary, "kernels.total.max"),
            "kernelTimeMean": get(summary, "kernels.total.mean"),
            "outsideKernelsMax": get(summary, "time.simulation.outsideKernels.max"),
            "outsideKernelsMean": get(summary, "time.simulation.outsideKernels.mean"),
            "loadImbalance": get(summary, "kernels.loadImbalance"),
            "cellUpdates": get(summary, "simulation.cellUpdates"),
            "duplicateCellUpdates": get(summary, "simulation.duplicateCellUpdates"),
            "cellUpdatesPerSecond": get(summary, "performance.cellUpdatesPerSecond"),
            "hardwareFlopPerSecond": get(summary, "performance.hardwareFlopPerSecond"),
            "nonzeroFlopPerSecond": get(summary, "performance.nonzeroFlopPerSecond"),
            "hostMemoryPeak": get(summary, "memory.atEnd.hostPeak.max"),
            "version": get(summary, "build.version"),
            "commit": get(summary, "build.commit"),
        }
    )
    if row["cells"] and run.ranks:
        row["cellsPerRank"] = row["cells"] / run.ranks
    ranks_reported = get(summary, "parallel.ranks")
    if ranks_reported is not None and ranks_reported != run.ranks:
        problems.append(f"it ran on {ranks_reported} ranks, not {run.ranks}")
    if campaign.machine.scheduler != "local" and row["hosts"] not in (None, run.nodes):
        problems.append(f"it ran on {row['hosts']} nodes, not {run.nodes}")

    def tagged(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
        kept = []
        for entry in rows:
            if entry.get("runId") not in (None, attempt["runId"]):
                continue
            kept.append(
                dict(
                    entry,
                    run=run.name,
                    build=run.build,
                    variant=run.variant,
                    repetition=run.repetition,
                    nodes=run.nodes,
                )
            )
        return kept

    epochs = tagged(read_csv(output / "data-epochs.csv"))
    steady = steady_state(epochs, campaign.job.warmup_epochs)
    row["steadyEpochs"], row["steadyWallTime"], row["steadySimulatedTime"] = steady
    ranks = tagged(read_csv(output / "data-ranks.csv"))
    pinning = pinning_problems(
        [
            r
            for r in read_csv(output / "data-threadPinning.csv")
            if r.get("runId") == attempt["runId"]
        ],
        run.layout.threads_per_rank,
    )
    row["pinningProblems"] = len(pinning)
    problems += pinning[:5]
    if problems:
        row["problems"] = "; ".join(problems)
    if status != COMPLETED:
        return row, [], []
    return row, ranks, epochs


def collect(
    campaign: Campaign, plan: Plan, directory: Directory, state: State
) -> dict[str, Path]:
    """Gathers the runs of @p plan into the tables of results/; returns their paths."""
    jobs = {job_id: info.get("state") for job_id, info in state.jobs.items()}
    rows, rank_rows, epoch_rows = [], [], []
    for run in plan.runs:
        attempt, status, count = chosen_attempt(run, directory, state, jobs)
        row, ranks, epochs = run_row(campaign, run, attempt, status, count, directory)
        rows.append(row)
        rank_rows += ranks
        epoch_rows += epochs
    identifying = ["run", "build", "variant", "repetition", "nodes"]

    def columns(entries: list[dict[str, Any]]) -> list[str]:
        names = list(identifying)
        for entry in entries:
            for key in entry:
                if key not in names:
                    names.append(key)
        return names

    results = directory.results
    return {
        "runs": write_csv(results / "runs.csv", RUN_COLUMNS, rows),
        "ranks": write_csv(results / "ranks.csv", columns(rank_rows), rank_rows),
        "epochs": write_csv(results / "epochs.csv", columns(epoch_rows), epoch_rows),
    }
