# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The job scripts of a campaign.

A job runs the attempts at the runs of one node count, one after the other.
Each attempt runs in a subshell of its own, with the environment of its
variant, and leaves in its directory what the campaign reads back: when it
started (started), how its check went (check.json, check.log), its log
(seissol.log), how it ended (exitcode, and check-failed if its check failed),
and the output of SeisSol (output/).
"""

from __future__ import annotations

import re
import shlex
import sys
from pathlib import Path
from typing import Any

from .campaign import Layout as Directory
from .campaign import utc_now
from .config import Campaign, Variant, format_duration
from .plan import Run
from .scheduler import JobSpec, Scheduler

PLACEHOLDERS = (
    "nodes",
    "ranks",
    "ranks_per_node",
    "cpus_per_rank",
    "threads_per_rank",
    "gpus_per_rank",
)
TOOL = Path(__file__).resolve().parent.parent / "scaling.py"


def fill_launcher(template: str, run: Run) -> str:
    """The launcher of @p run: the template with the values of the run put in.

    Only the placeholders of the run are replaced, so that the shell variables a
    launcher may use, e.g. ${CPU_BIND}, stay as they are.
    """
    values = {
        "nodes": run.nodes,
        "ranks": run.ranks,
        "ranks_per_node": run.layout.ranks_per_node,
        "cpus_per_rank": run.layout.cpus_per_rank,
        "threads_per_rank": run.layout.threads_per_rank,
        "gpus_per_rank": run.layout.gpus_per_rank,
    }
    pattern = re.compile(r"(?<!\$)\{(" + "|".join(PLACEHOLDERS) + r")\}")
    return pattern.sub(lambda match: str(values[match.group(1)]), template)


def run_overrides(campaign: Campaign, variant: Variant, prefix: str) -> dict[str, str]:
    """The --set overrides of a run of @p variant, which writes its output at @p prefix."""
    overrides = dict(campaign.case.overrides)
    if campaign.case.epoch_interval is not None:
        overrides["output.epochinterval"] = repr(campaign.case.epoch_interval)
    overrides.update(variant.overrides)
    overrides["output.outputfile"] = prefix
    return overrides


def seissol_arguments(campaign: Campaign, overrides: dict[str, str]) -> list[str]:
    return [str(campaign.case.parameter_file)] + [
        f"--set={name}={value}" for name, value in overrides.items()
    ]


def export_lines(env: dict[str, str]) -> list[str]:
    return [f"export {name}={shlex.quote(value)}" for name, value in env.items()]


def with_time_limit(launcher: str, scheduler: str, seconds: int | None) -> str:
    """The launcher, limited to @p seconds where the launcher can be."""
    if seconds is None:
        return launcher
    if scheduler == "slurm" and launcher.split()[:1] == ["srun"]:
        return launcher.replace("srun", f"srun --time={format_duration(seconds)}", 1)
    return f"timeout {seconds} {launcher}"


def run_block(
    campaign: Campaign,
    directory: Directory,
    run: Run,
    attempt: dict[str, Any],
    binary: Path,
    index: int,
    count: int,
) -> list[str]:
    variant = next(v for v in campaign.variants if v.name == run.variant)
    attempt_dir = directory.attempt(run.name, attempt["attempt"])
    overrides = run_overrides(campaign, variant, str(attempt_dir / "output" / "data"))
    arguments = seissol_arguments(campaign, overrides)
    env = {"OMP_NUM_THREADS": str(run.layout.threads_per_rank)}
    env.update(variant.env)
    env["SEISSOL_RUN_ID"] = attempt["runId"]
    launcher = fill_launcher(campaign.machine.launcher, run)
    if variant.launcher_args:
        launcher += " " + " ".join(variant.launcher_args)
    launcher = with_time_limit(
        launcher, campaign.machine.scheduler, campaign.job.run_time_limit
    )
    wrapper = (
        f" {shlex.quote(campaign.machine.wrapper)}" if campaign.machine.wrapper else ""
    )
    lines = [
        "",
        f"# run {index} of {count}: {run.name}, attempt {attempt['attempt']}",
        "(",
        f"  attempt={shlex.quote(str(attempt_dir))}",
        f"  seissol={shlex.quote(str(binary))}",
        "  arguments=(",
    ]
    lines += [f"    {shlex.quote(argument)}" for argument in arguments]
    lines += [
        "  )",
        '  mkdir -p "$attempt/output"',
        f"  cd {shlex.quote(str(campaign.case.parameter_file.parent))} || exit 1",
    ]
    lines += ["  " + line for line in export_lines(env)]
    lines += [
        '  { date -u +%Y-%m-%dT%H:%M:%SZ; hostname; } > "$attempt/started"',
        f"  echo \"run {run.name}, attempt {attempt['attempt']}, $(date -u +%H:%M:%S)\"",
    ]
    command = f'{launcher}{wrapper} "$seissol" "${{arguments[@]}}"'
    if campaign.job.check_in_job:
        lines += [
            '  if "$seissol" --check="$attempt/check.json" "${arguments[@]}" '
            '> "$attempt/check.log" 2>&1; then',
            f'    {command} > "$attempt/seissol.log" 2>&1',
            '    echo $? > "$attempt/exitcode"',
            "  else",
            '    echo $? > "$attempt/check-failed"',
            '    echo 1 > "$attempt/exitcode"',
            "  fi",
        ]
    else:
        lines += [
            f'  {command} > "$attempt/seissol.log" 2>&1',
            '  echo $? > "$attempt/exitcode"',
        ]
    lines += [
        '  echo "  exit code $(cat "$attempt/exitcode"), $(date -u +%H:%M:%S)"',
        ")",
    ]
    return lines


def job_spec(
    campaign: Campaign, directory: Directory, name: str, runs: list[Run]
) -> JobSpec:
    # the layout of one of the runs, which fits a node, where the largest of each would not;
    # every run states its own to the launcher
    widest = max(
        runs, key=lambda run: (run.layout.ranks_per_node, run.layout.cpus_per_rank)
    ).layout
    gpus = max(run.layout.ranks_per_node * run.layout.gpus_per_rank for run in runs)
    return JobSpec(
        name=f"{campaign.name}-{name}",
        nodes=max(run.nodes for run in runs),
        time_limit=campaign.job.time_limit,
        ranks_per_node=widest.ranks_per_node,
        cpus_per_rank=widest.cpus_per_rank,
        gpus_per_node=gpus,
        log=directory.jobs / f"{name}-%j.out",
        workdir=directory.root,
    )


def measurement_script(
    campaign: Campaign,
    directory: Directory,
    scheduler: Scheduler,
    name: str,
    entries: list[tuple[Run, dict[str, Any]]],
    binaries: dict[str, Path],
) -> str:
    """The job script that runs the attempts @p entries, of the job @p name."""
    runs = [run for run, _ in entries]
    spec = job_spec(campaign, directory, name, runs)
    machine = campaign.machine
    lines = ["#!/bin/bash"]
    lines += scheduler.header(spec)
    lines += [
        "#",
        f"# Strong scaling campaign {campaign.name}, job {name}: {len(entries)} run(s) "
        f"on {spec.nodes} node(s).",
        f"# Written by {TOOL.name} at {utc_now()}; every run writes into its directory",
        f"# under {directory.runs}.",
        "",
        "set -u",
    ]
    if machine.setup:
        lines += ["", "# the setup of the machine"] + machine.setup
    if machine.env:
        lines += ["", "# the environment of the machine"] + export_lines(machine.env)
    if campaign.case.env:
        lines += ["", "# the environment of the case"] + export_lines(campaign.case.env)
    lines += [
        "",
        'echo "job ${SLURM_JOB_ID:-local} on $(hostname), $(date -u +%Y-%m-%dT%H:%M:%SZ)"',
    ]
    for index, (run, attempt) in enumerate(entries, start=1):
        lines += run_block(
            campaign, directory, run, attempt, binaries[run.build], index, len(entries)
        )
    lines += ["", 'echo "job done, $(date -u +%Y-%m-%dT%H:%M:%SZ)"', ""]
    return "\n".join(lines)


def report_script(
    campaign: Campaign, directory: Directory, scheduler: Scheduler
) -> str:
    """The job script that writes the report once the runs have ended."""
    spec = JobSpec(
        name=f"{campaign.name}-report",
        nodes=1,
        time_limit=30 * 60,
        ranks_per_node=1,
        cpus_per_rank=1,
        gpus_per_node=0,
        log=directory.jobs / "report-%j.out",
        workdir=directory.root,
        exclusive=False,
    )
    python = campaign.machine.python or sys.executable
    lines = ["#!/bin/bash"]
    lines += scheduler.header(spec)
    lines += ["", "set -u"]
    if campaign.machine.setup:
        lines += campaign.machine.setup
    lines += [
        f"cd {shlex.quote(str(directory.root))}",
        f"{python} {shlex.quote(str(TOOL))} report {shlex.quote(str(directory.root))}",
        "",
    ]
    return "\n".join(lines)
