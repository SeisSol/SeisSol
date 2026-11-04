# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Checks of a campaign before its jobs are submitted.

What can be found out before waiting in a queue: whether the layout of the
ranks fits a node, whether the scheduler and the launcher are there, and, by
`SeisSol --check` of every build and variant, whether the parameters, the files
they name, the mesh and the output directory are fine. The check of SeisSol
also gives the size of the mesh, which the scaling guards need.
"""

from __future__ import annotations

import os
import signal
import subprocess  # nosec B404
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .campaign import Layout as Directory
from .campaign import read_json
from .config import Campaign, ConfigError, partition_for
from .jobscript import run_overrides, seissol_arguments
from .scheduler import Scheduler


@dataclass
class Finding:
    level: str  # error, warning or note
    text: str


@dataclass
class CheckResult:
    findings: list[Finding] = field(default_factory=list)
    cells: int | None = None
    reports: dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not any(f.level == "error" for f in self.findings)

    def add(self, level: str, text: str) -> None:
        self.findings.append(Finding(level, text))


def check_layouts(campaign: Campaign, result: CheckResult) -> None:
    machine = campaign.machine
    cpus_per_node = machine.cores_per_node * machine.threads_per_core
    for variant in campaign.variants:
        layout = variant.layout(machine)
        name = f"variant {variant.name}"
        taken = layout.ranks_per_node * layout.cpus_per_rank
        if taken > cpus_per_node:
            result.add(
                "error",
                f"{name}: {layout.ranks_per_node} ranks of {layout.cpus_per_rank} CPUs take "
                f"{taken} CPUs of a node, which has {cpus_per_node}",
            )
        if layout.threads_per_rank > layout.cpus_per_rank:
            result.add(
                "error",
                f"{name}: {layout.threads_per_rank} threads per rank do not fit "
                f"{layout.cpus_per_rank} CPUs per rank",
            )
        commthread = variant.env.get(
            "SEISSOL_COMMTHREAD", machine.env.get("SEISSOL_COMMTHREAD", "1")
        )
        if layout.threads_per_rank >= layout.cpus_per_rank and commthread != "0":
            result.add(
                "warning",
                f"{name}: the {layout.threads_per_rank} threads of a rank take all its CPUs, "
                "and none is left for the communication thread; set SEISSOL_COMMTHREAD=0 "
                "or give a rank a CPU more",
            )
        gpus = layout.ranks_per_node * layout.gpus_per_rank
        if gpus > machine.gpus_per_node:
            result.add(
                "error",
                f"{name}: {layout.ranks_per_node} ranks of {layout.gpus_per_rank} GPUs need "
                f"{gpus} GPUs, but a node has {machine.gpus_per_node}",
            )
    wrapper = machine.wrapper
    if wrapper and os.path.isabs(wrapper):
        if not os.path.exists(wrapper):
            result.add("error", f"the wrapper {wrapper} of the machine does not exist")
        elif not os.access(wrapper, os.X_OK):
            result.add(
                "error", f"the wrapper {wrapper} of the machine is not executable"
            )


def check_partitions(
    campaign: Campaign, node_counts: list[int], result: CheckResult
) -> None:
    for nodes in node_counts:
        try:
            partition_for(campaign.machine.slurm, nodes)
        except ConfigError as error:
            result.add("error", str(error))


def check_tools(scheduler: Scheduler, result: CheckResult, submitting: bool) -> None:
    missing = scheduler.missing_tools()
    if missing:
        result.add(
            "error" if submitting else "warning",
            f"the {scheduler.name} commands {', '.join(missing)} are not available here",
        )


def _last_error(log: Path) -> str:
    try:
        lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    errors = [line for line in lines if " error " in line]
    return errors[-1].split(" : ", 1)[-1].strip() if errors else ""


def seissol_check(
    campaign: Campaign,
    directory: Directory,
    build: str,
    binary: Path,
    variant_name: str,
    result: CheckResult,
    timeout: int = 600,
) -> dict[str, Any] | None:
    """Runs SeisSol --check for a run of @p build and @p variant_name; returns its report."""
    variant = next(v for v in campaign.variants if v.name == variant_name)
    key = f"{build}-{variant_name}"
    directory.checks.mkdir(parents=True, exist_ok=True)
    report_path = directory.checks / f"{key}.json"
    log_path = directory.checks / f"{key}.log"
    report_path.unlink(missing_ok=True)
    overrides = run_overrides(campaign, variant, str(directory.checks / key / "data"))
    command = [str(binary), f"--check={report_path}"] + seissol_arguments(
        campaign, overrides
    )
    env = dict(os.environ)
    env.update(campaign.machine.env)
    env.update(campaign.case.env)
    env.update(variant.env)
    where = f"the check of {build}/{variant_name}"
    try:
        with log_path.open("w", encoding="utf-8") as stream:
            completed = subprocess.run(  # nosec B603
                command,
                cwd=campaign.case.parameter_file.parent,
                env=env,
                stdout=stream,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
    except subprocess.TimeoutExpired:
        result.add("warning", f"{where} did not finish in {timeout} s; see {log_path}")
        return None
    except OSError as error:
        result.add("error", f"{where}: cannot run {binary}: {error.strerror}")
        return None
    if completed.returncode < 0:
        number = -completed.returncode
        try:
            name = signal.Signals(number).name
        except ValueError:
            name = f"signal {number}"
        result.add(
            "warning",
            f"{where}: {binary.name} does not run on this node ({name}); built for the "
            "compute nodes, it may not run on a login node. The jobs check every run "
            "before running it.",
        )
        return None
    if not report_path.exists():
        result.add(
            "error",
            f"{where} failed with exit code {completed.returncode} and wrote no report; "
            f"see {log_path}",
        )
        return None
    report = read_json(report_path)
    if report.get("status") != "passed":
        for problem in report.get("problems", []):
            result.add("error", f"{where}: {problem}")
        if report.get("stoppedAt"):
            message = _last_error(log_path)
            result.add(
                "error",
                f"{where} stopped at {report['stoppedAt']}"
                + (f": {message}" if message else "")
                + f" (see {log_path})",
            )
        if not report.get("problems") and not report.get("stoppedAt"):
            result.add("error", f"{where} failed; see {log_path}")
    for unknown in report.get("unknownFields", []):
        result.add(
            "warning",
            f"{where}: the parameter file has {unknown}, which SeisSol does not read",
        )
    enabled = report.get("output", {}).get("enabled", [])
    if enabled:
        result.add(
            "note",
            f"{where}: the runs write {', '.join(enabled)}, which their times include",
        )
    return report


def check_campaign(
    campaign: Campaign,
    directory: Directory,
    scheduler: Scheduler,
    binaries: dict[str, Path],
    submitting: bool = False,
    run_seissol: bool = True,
    node_counts: list[int] | None = None,
) -> CheckResult:
    """All checks of @p campaign, with the executables @p binaries of its builds."""
    result = CheckResult()
    if not campaign.case.parameter_file.exists():
        result.add(
            "error", f"the parameter file {campaign.case.parameter_file} does not exist"
        )
    check_layouts(campaign, result)
    if campaign.machine.scheduler == "slurm":
        check_partitions(campaign, node_counts or [], result)
    check_tools(scheduler, result, submitting)
    for build in campaign.builds:
        if build.name not in binaries:
            result.add(
                "note", f"build {build.name} is not built yet; 'build' builds it"
            )
    if not run_seissol or not campaign.case.parameter_file.exists():
        return result
    cells = set()
    for build in campaign.builds:
        binary = binaries.get(build.name)
        if binary is None:
            continue
        for variant in campaign.variants:
            report = seissol_check(
                campaign, directory, build.name, binary, variant.name, result
            )
            if report is None:
                continue
            result.reports[f"{build.name}/{variant.name}"] = report
            mesh_cells = report.get("mesh", {}).get("cells")
            if mesh_cells:
                cells.add(int(mesh_cells))
    if len(cells) > 1:
        result.add(
            "warning",
            f"the checks found meshes of different sizes: {', '.join(map(str, sorted(cells)))} "
            "cells; the variants should run the same mesh",
        )
    if cells:
        result.cells = max(cells)
    return result
