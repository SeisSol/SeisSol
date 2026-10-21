# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The batch systems jobs are submitted to: Slurm, and the local machine.

A scheduler writes the header of a job script, submits it with dependencies
on other jobs, and tells the state of jobs as pending, running or ended.
Another batch system is another subclass of Scheduler.
"""

from __future__ import annotations

import itertools
import os
import shutil
import subprocess  # nosec B404
from dataclasses import dataclass
from pathlib import Path

from .config import Machine, format_duration, partition_for


class SchedulerError(Exception):
    """A job could not be submitted, or its state not be found out."""


@dataclass
class JobSpec:
    name: str
    nodes: int
    #: seconds
    time_limit: int
    ranks_per_node: int
    cpus_per_rank: int
    gpus_per_node: int
    #: where the output of the job goes
    log: Path
    workdir: Path
    #: whether the job takes its nodes for itself
    exclusive: bool = True


class Scheduler:
    name = ""

    def header(self, job: JobSpec) -> list[str]:
        raise NotImplementedError

    def submit(self, script: Path, after: list[str], kind: str = "afterok") -> str:
        """Submits @p script to start once the jobs @p after are done; returns its ID.

        kind is afterok (once they succeeded) or afterany (once they ended).
        """
        raise NotImplementedError

    def states(self, job_ids: list[str]) -> dict[str, str]:
        """pending, running or ended per job; jobs that are not known are left out."""
        raise NotImplementedError

    def cancel(self, job_ids: list[str]) -> None:
        raise NotImplementedError

    def missing_tools(self) -> list[str]:
        """The commands the scheduler needs that are not there."""
        return []


SLURM_PENDING = {
    "PENDING",
    "CONFIGURING",
    "REQUEUED",
    "REQUEUE_HOLD",
    "REQUEUE_FED",
    "RESV_DEL_HOLD",
}
SLURM_RUNNING = {
    "RUNNING",
    "COMPLETING",
    "STAGE_OUT",
    "SIGNALING",
    "RESIZING",
    "SUSPENDED",
}


def slurm_state(text: str) -> str:
    """pending, running or ended of a state as squeue and sacct write it."""
    word = text.strip().split()[0].rstrip("+").upper() if text.strip() else ""
    if word in SLURM_PENDING:
        return "pending"
    if word in SLURM_RUNNING:
        return "running"
    return "ended"


class SlurmScheduler(Scheduler):
    name = "slurm"

    def __init__(self, machine: Machine):
        self.options = machine.slurm

    def header(self, job: JobSpec) -> list[str]:
        lines = [
            f"#SBATCH --job-name={job.name}",
            f"#SBATCH --nodes={job.nodes}",
            f"#SBATCH --ntasks-per-node={job.ranks_per_node}",
        ]
        if self.options.get("cpus_per_task", True):
            lines.append(f"#SBATCH --cpus-per-task={job.cpus_per_rank}")
        lines += [
            f"#SBATCH --time={format_duration(job.time_limit)}",
            f"#SBATCH --output={job.log}",
            f"#SBATCH --chdir={job.workdir}",
            # an attempt that fails is resubmitted by the campaign, as a new attempt
            "#SBATCH --no-requeue",
        ]
        partition = partition_for(self.options, job.nodes)
        if partition:
            lines.append(f"#SBATCH --partition={partition}")
        for key in ("account", "qos", "constraint"):
            value = self.options.get(key)
            if value:
                lines.append(f"#SBATCH --{key}={value}")
        if job.exclusive and self.options.get("exclusive", True):
            lines.append("#SBATCH --exclusive")
        if job.gpus_per_node:
            lines.append(f"#SBATCH --gpus-per-node={job.gpus_per_node}")
        for option in self.options.get("extra", []):
            option = option.strip()
            if option.startswith("#SBATCH"):
                lines.append(option)
            else:
                lines.append(f"#SBATCH {option}")
        return lines

    def missing_tools(self) -> list[str]:
        return [
            tool
            for tool in ("sbatch", "squeue", "scancel")
            if shutil.which(tool) is None
        ]

    def submit(self, script: Path, after: list[str], kind: str = "afterok") -> str:
        command = ["sbatch", "--parsable"]
        if after:
            command.append(f"--dependency={kind}:{':'.join(after)}")
        command.append(str(script))
        try:
            result = subprocess.run(  # nosec B603
                command, check=False, capture_output=True, text=True
            )
        except OSError as error:
            raise SchedulerError(f"cannot run sbatch: {error.strerror}") from error
        if result.returncode != 0:
            raise SchedulerError(
                f"sbatch refused {script}: {(result.stderr or result.stdout).strip()}"
            )
        job_id = result.stdout.strip().split(";")[0]
        if not job_id:
            raise SchedulerError(f"sbatch gave no job ID for {script}")
        return job_id

    def states(self, job_ids: list[str]) -> dict[str, str]:
        if not job_ids:
            return {}
        states: dict[str, str] = {}
        wanted = set(job_ids)
        user = os.environ.get("USER") or os.environ.get("LOGNAME")
        command = ["squeue", "--noheader", "--format=%i|%T"]
        if user:
            command.append(f"--user={user}")
        try:
            result = subprocess.run(  # nosec B603
                command, check=False, capture_output=True, text=True
            )
        except OSError as error:
            raise SchedulerError(f"cannot run squeue: {error.strerror}") from error
        if result.returncode != 0:
            # a queue that cannot be read says nothing about the jobs in it
            raise SchedulerError(
                f"squeue failed: {(result.stderr or result.stdout).strip()}"
            )
        for line in result.stdout.splitlines():
            job_id, _, state = line.partition("|")
            if job_id.strip() in wanted:
                states[job_id.strip()] = slurm_state(state)
        rest = [job_id for job_id in job_ids if job_id not in states]
        if rest and shutil.which("sacct") is not None:
            result = subprocess.run(  # nosec B603
                [
                    "sacct",
                    "--noheader",
                    "--parsable2",
                    "--allocations",
                    "--format=JobID,State",
                    f"--jobs={','.join(rest)}",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                for line in result.stdout.splitlines():
                    job_id, _, state = line.partition("|")
                    if job_id.strip() in wanted:
                        states[job_id.strip()] = slurm_state(state)
        # a job that is neither queued nor accounted for has ended long ago
        for job_id in rest:
            states.setdefault(job_id, "ended")
        return states

    def cancel(self, job_ids: list[str]) -> None:
        if job_ids:
            subprocess.run(["scancel", *job_ids], check=False)  # nosec B603 B607


class LocalScheduler(Scheduler):
    """Runs a job on this machine, at once, when it is submitted."""

    name = "local"
    _counter = itertools.count(1)

    def __init__(self, machine: Machine):
        self.machine = machine

    def header(self, job: JobSpec) -> list[str]:
        return [
            f"# job {job.name}: {job.nodes} node(s) of {job.ranks_per_node} rank(s), "
            f"time limit {format_duration(job.time_limit)}, run on this machine"
        ]

    def missing_tools(self) -> list[str]:
        launcher = (
            self.machine.launcher.split()[0] if self.machine.launcher else "mpirun"
        )
        return [] if shutil.which(launcher) else [launcher]

    def submit(self, script: Path, after: list[str], kind: str = "afterok") -> str:
        job_id = f"local-{os.getpid()}-{next(self._counter)}"
        log = script.with_suffix(".log")
        with log.open("w", encoding="utf-8") as stream:
            subprocess.run(  # nosec B603 B607
                ["bash", str(script)],
                check=False,
                stdout=stream,
                stderr=subprocess.STDOUT,
            )
        return job_id

    def states(self, job_ids: list[str]) -> dict[str, str]:
        return {job_id: "ended" for job_id in job_ids}

    def cancel(self, job_ids: list[str]) -> None:
        return None


def scheduler_for(machine: Machine) -> Scheduler:
    if machine.scheduler == "slurm":
        return SlurmScheduler(machine)
    return LocalScheduler(machine)
