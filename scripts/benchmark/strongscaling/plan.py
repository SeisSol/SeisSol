# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The runs of a campaign: which node counts, and in which order.

Every node count gets a job of its own, which runs all builds and variants at
that node count. Within the job, the repetitions take turns (interleaved): the
first repetition of every build and variant, then the second, each time in a
rotated order, so that a slow phase of the machine does not fall on one
variant alone.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from .config import Campaign, ConfigError, Layout

PLAN_SCHEMA = "seissol-scaling-plan/1"


@dataclass
class Run:
    name: str
    nodes: int
    build: str
    variant: str
    repetition: int
    layout: Layout
    job: str
    #: a hash of all that the run depends on; a run whose definition changed is
    #: measured anew
    definition: str

    @property
    def ranks(self) -> int:
        return self.nodes * self.layout.ranks_per_node

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "nodes": self.nodes,
            "ranks": self.ranks,
            "build": self.build,
            "variant": self.variant,
            "repetition": self.repetition,
            "layout": self.layout.to_dict(),
            "job": self.job,
            "definition": self.definition,
        }


@dataclass
class JobPlan:
    name: str
    nodes: int
    runs: list[str]


@dataclass
class Plan:
    campaign: str
    cells: int | None
    node_counts: list[int]
    runs: list[Run]
    jobs: list[JobPlan]
    notes: list[str] = field(default_factory=list)

    def run(self, name: str) -> Run:
        for run in self.runs:
            if run.name == name:
                return run
        raise KeyError(name)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": PLAN_SCHEMA,
            "campaign": self.campaign,
            "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "cells": self.cells,
            "nodeCounts": self.node_counts,
            "notes": self.notes,
            "runs": [run.to_dict() for run in self.runs],
            "jobs": [
                {"name": job.name, "nodes": job.nodes, "runs": job.runs}
                for job in self.jobs
            ],
        }


def candidate_node_counts(
    campaign: Campaign, cells: int | None
) -> tuple[list[int], list[str]]:
    """The node counts of the scaling settings, before the guards."""
    scaling = campaign.scaling
    if scaling.nodes is not None:
        return list(scaling.nodes), []
    low, high, factor = scaling.node_range
    notes = []
    if low is None:
        if cells is None or scaling.max_cells_per_node is None:
            raise ConfigError(
                "scaling.nodes.min: 'auto' takes the smallest node count the mesh fits on, "
                "which needs the cells of the mesh (from 'check', or case.cells) and "
                "scaling.max_cells_per_node"
            )
        low = max(1, math.ceil(cells / scaling.max_cells_per_node))
        notes.append(
            f"the smallest node count is {low}: {cells} cells at most "
            f"{scaling.max_cells_per_node:g} per node"
        )
    counts = []
    nodes = low
    while nodes <= high:
        counts.append(nodes)
        nodes *= factor
    if not counts:
        raise ConfigError(f"scaling.nodes: no node count between {low} and {high}")
    return counts, notes


def guarded(
    campaign: Campaign, nodes: int, layout: Layout, cells: int | None
) -> str | None:
    """Why a run on @p nodes nodes is left out, if it is."""
    scaling = campaign.scaling
    if cells is None:
        return None
    if (
        scaling.max_cells_per_node is not None
        and cells / nodes > scaling.max_cells_per_node
    ):
        return (
            f"{cells / nodes:.0f} cells per node exceed "
            f"scaling.max_cells_per_node = {scaling.max_cells_per_node:g}"
        )
    ranks = nodes * layout.ranks_per_node
    if (
        scaling.min_cells_per_rank is not None
        and cells / ranks < scaling.min_cells_per_rank
    ):
        return (
            f"{cells / ranks:.0f} cells per rank fall below "
            f"scaling.min_cells_per_rank = {scaling.min_cells_per_rank:g}"
        )
    return None


def definition_hash(
    campaign: Campaign,
    nodes: int,
    build: str,
    variant: str,
    fingerprint: str | None = None,
) -> str:
    """A hash of everything a run depends on, the executable among it (@p fingerprint)."""
    build_spec = next(b for b in campaign.builds if b.name == build)
    variant_spec = next(v for v in campaign.variants if v.name == variant)
    machine = campaign.machine
    payload = {
        "nodes": nodes,
        "case": {
            "parameterFile": str(campaign.case.parameter_file),
            "set": campaign.case.overrides,
            "epochInterval": campaign.case.epoch_interval,
            "env": campaign.case.env,
        },
        "build": {
            "binary": str(build_spec.binary) if build_spec.binary else None,
            "variants": (
                str(build_spec.variants_file) if build_spec.variants_file else None
            ),
            "preset": build_spec.preset,
            "toolchain": build_spec.toolchain,
            "set": build_spec.overrides,
            "select": build_spec.select,
            "fingerprint": fingerprint,
        },
        "variant": {
            "env": variant_spec.env,
            "set": variant_spec.overrides,
            "layout": variant_spec.layout(machine).to_dict(),
            "launcherArgs": variant_spec.launcher_args,
        },
        "machine": {
            "launcher": machine.launcher,
            "wrapper": machine.wrapper,
            "env": machine.env,
            "setup": machine.setup,
        },
    }
    text = json.dumps(payload, sort_keys=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def job_name(nodes: int) -> str:
    return f"n{nodes:04d}"


def run_name(nodes: int, build: str, variant: str, repetition: int) -> str:
    return f"{job_name(nodes)}-{build}-{variant}-r{repetition}"


def make_plan(
    campaign: Campaign,
    cells: int | None = None,
    fingerprints: dict[str, str] | None = None,
) -> Plan:
    """The runs and jobs of @p campaign, for a mesh of @p cells cells if known.

    @p fingerprints are those of the executables of the builds: a run of another
    executable than the one an attempt ran is another run.
    """
    fingerprints = fingerprints or {}
    cells = cells if cells is not None else campaign.case.cells
    counts, notes = candidate_node_counts(campaign, cells)
    if cells is None and (
        campaign.scaling.min_cells_per_rank is not None
        or campaign.scaling.max_cells_per_node is not None
    ):
        notes.append(
            "the cells of the mesh are unknown, so the guards of scaling are not "
            "applied; 'check' finds them"
        )
    combinations = [
        (build.name, variant)
        for build in campaign.builds
        for variant in campaign.variants
    ]
    runs: list[Run] = []
    jobs: list[JobPlan] = []
    kept_counts = []
    for nodes in counts:
        allowed = []
        for build, variant in combinations:
            layout = variant.layout(campaign.machine)
            reason = guarded(campaign, nodes, layout, cells)
            if reason is not None:
                notes.append(
                    f"{nodes} nodes, {build}/{variant.name}: left out, {reason}"
                )
                continue
            allowed.append((build, variant, layout))
        if not allowed:
            continue
        kept_counts.append(nodes)
        order = []
        repetitions = range(1, campaign.scaling.repetitions + 1)
        if campaign.scaling.order == "interleaved":
            for repetition in repetitions:
                shift = (repetition - 1) % len(allowed)
                rotated = allowed[shift:] + allowed[:shift]
                order += [(entry, repetition) for entry in rotated]
        else:
            order = [
                (entry, repetition) for entry in allowed for repetition in repetitions
            ]
        job = JobPlan(name=job_name(nodes), nodes=nodes, runs=[])
        for (build, variant, layout), repetition in order:
            run = Run(
                name=run_name(nodes, build, variant.name, repetition),
                nodes=nodes,
                build=build,
                variant=variant.name,
                repetition=repetition,
                layout=layout,
                job=job.name,
                definition=definition_hash(
                    campaign, nodes, build, variant.name, fingerprints.get(build)
                ),
            )
            runs.append(run)
            job.runs.append(run.name)
        jobs.append(job)
    if not runs:
        raise ConfigError(
            "the guards of scaling leave no run; see the notes: " + "; ".join(notes)
        )
    return Plan(
        campaign=campaign.name,
        cells=cells,
        node_counts=kept_counts,
        runs=runs,
        jobs=jobs,
        notes=notes,
    )


def describe(plan: Plan) -> list[str]:
    """The plan as lines of a table, for the terminal."""
    lines = []
    cells = f"{plan.cells} cells" if plan.cells is not None else "cells unknown"
    lines.append(
        f"campaign {plan.campaign}: {len(plan.runs)} runs in {len(plan.jobs)} jobs ({cells})"
    )
    lines.append(f"{'nodes':>6} {'ranks':>6} {'cells/rank':>11}  runs")
    for job in plan.jobs:
        job_runs = [plan.run(name) for name in job.runs]
        ranks = sorted({run.ranks for run in job_runs})
        per_rank = (
            ", ".join(f"{plan.cells / r:.0f}" for r in ranks)
            if plan.cells is not None
            else "-"
        )
        lines.append(
            f"{job.nodes:>6} {'/'.join(map(str, ranks)):>6} {per_rank:>11}  {len(job_runs)}"
        )
    for note in plan.notes:
        lines.append(f"note: {note}")
    return lines
