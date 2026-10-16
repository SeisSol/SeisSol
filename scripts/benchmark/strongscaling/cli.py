# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The commands of scaling.py."""

from __future__ import annotations

import argparse
import contextlib
import sys
from pathlib import Path
from typing import Any

from . import plan as planning
from .build import BuildError, build_all, known_binaries, known_fingerprints
from .campaign import (
    ACTIVE_STATES,
    COMPLETED,
    Layout,
    State,
    attempt_status,
    locked,
    utc_now,
    write_json,
    write_text_atomically,
)
from .check import CheckResult, check_campaign
from .config import Campaign, ConfigError, load_campaign
from .jobscript import measurement_script
from .scheduler import Scheduler, SchedulerError, scheduler_for


class CommandError(Exception):
    """A command that cannot go on; its message says why."""


def log(message: str) -> None:
    print(message, flush=True)


def open_campaign(args: argparse.Namespace) -> tuple[Campaign, Layout, State]:
    """The campaign the command names, held for the command alone."""
    campaign = load_campaign(Path(args.campaign))
    directory = Layout(campaign.directory)
    args.stack.enter_context(locked(directory.root))
    write_json(
        directory.config_copy,
        {
            "path": str(campaign.path),
            "config": campaign.raw,
            "machine": campaign.machine.to_dict(),
        },
    )
    return campaign, directory, State(directory.state)


def known_cells(
    campaign: Campaign, state: State, given: int | None = None
) -> int | None:
    if given is not None:
        return given
    if campaign.case.cells is not None:
        return campaign.case.cells
    cells = [
        entry.get("cells") for entry in state.checks.values() if entry.get("cells")
    ]
    return max(cells) if cells else None


def make_plan(
    campaign: Campaign, state: State, cells: int | None = None
) -> planning.Plan:
    return planning.make_plan(
        campaign,
        known_cells(campaign, state, cells),
        known_fingerprints(campaign, state.builds),
    )


def write_plan(
    campaign: Campaign, directory: Layout, state: State, cells: int | None = None
) -> planning.Plan:
    plan = make_plan(campaign, state, cells)
    write_json(directory.plan, plan.to_dict())
    return plan


def print_findings(result: CheckResult) -> None:
    for finding in result.findings:
        log(f"{finding.level}: {finding.text}")


def record_checks(state: State, result: CheckResult) -> None:
    for key, report in result.reports.items():
        state.checks[key] = {
            "status": report.get("status"),
            "cells": report.get("mesh", {}).get("cells"),
            "checked": utc_now(),
        }


# ---------------------------------------------------------------------------
# plan, check, build
# ---------------------------------------------------------------------------


def cmd_plan(args: argparse.Namespace) -> int:
    campaign, directory, state = open_campaign(args)
    plan = write_plan(campaign, directory, state, args.cells)
    for line in planning.describe(plan):
        log(line)
    log(f"wrote {directory.plan}")
    return 0


def run_checks(
    campaign: Campaign,
    directory: Layout,
    state: State,
    binaries: dict[str, Path],
    submitting: bool,
    run_seissol: bool = True,
) -> CheckResult:
    scheduler = scheduler_for(campaign.machine)
    try:
        counts = make_plan(campaign, state).node_counts
    except ConfigError:
        counts = []
    result = check_campaign(
        campaign,
        directory,
        scheduler,
        binaries,
        submitting=submitting,
        run_seissol=run_seissol,
        node_counts=counts,
    )
    record_checks(state, result)
    state.save()
    return result


def cmd_check(args: argparse.Namespace) -> int:
    campaign, directory, state = open_campaign(args)
    binaries = known_binaries(campaign, state.builds)
    result = run_checks(
        campaign,
        directory,
        state,
        binaries,
        submitting=False,
        run_seissol=not args.no_seissol,
    )
    print_findings(result)
    plan = write_plan(campaign, directory, state)
    for line in planning.describe(plan):
        log(line)
    if result.ok:
        log("the check passed")
        return 0
    log("the check failed")
    return 1


def cmd_build(args: argparse.Namespace) -> int:
    campaign, directory, state = open_campaign(args)
    try:
        binaries = build_all(campaign, directory, state.builds, log=log)
    finally:
        state.save()
    for name, binary in binaries.items():
        log(f"{name}: {binary}")
    return 0


# ---------------------------------------------------------------------------
# submit, status, cancel
# ---------------------------------------------------------------------------


def job_states(scheduler: Scheduler, state: State) -> dict[str, str]:
    """The states of the jobs of the campaign, with those that ended recorded as such."""
    open_jobs = [
        job_id for job_id, job in state.jobs.items() if job.get("state") != "ended"
    ]
    try:
        states = scheduler.states(open_jobs)
    except SchedulerError as error:
        log(f"warning: {error}")
        states = {}
    for job_id, job_state in states.items():
        state.jobs[job_id]["state"] = job_state
    return {job_id: job.get("state") for job_id, job in state.jobs.items()}


def run_states(
    plan: planning.Plan, directory: Layout, state: State, jobs: dict[str, str]
) -> dict[str, tuple[dict[str, Any] | None, str]]:
    """The latest attempt at every run of the plan, and its state."""
    states = {}
    for run in plan.runs:
        attempts = [
            a for a in state.attempts(run.name) if a["definition"] == run.definition
        ]
        if not attempts:
            states[run.name] = (None, "planned")
            continue
        attempt = attempts[-1]
        job_state = jobs.get(attempt.get("job")) if attempt.get("job") else None
        states[run.name] = (
            attempt,
            attempt_status(directory.attempt(run.name, attempt["attempt"]), job_state),
        )
    return states


def open_runs(
    campaign: Campaign, plan: planning.Plan, state: State, states: dict[str, Any]
) -> tuple[list[planning.Run], list[str]]:
    """The runs to attempt (again), and why the others that did not complete are not."""
    runs = []
    notes = []
    for run in plan.runs:
        attempt, status = states[run.name]
        if status == COMPLETED or status in ACTIVE_STATES or status == "unknown":
            continue
        tried = len(
            [a for a in state.attempts(run.name) if a["definition"] == run.definition]
        )
        if status == "check-failed":
            notes.append(f"{run.name}: its check failed; it is not attempted again")
            continue
        if tried >= campaign.job.max_attempts:
            notes.append(
                f"{run.name}: {tried} attempts failed; it is not attempted again"
            )
            continue
        runs.append(run)
    return runs, notes


def stale_jobs(plan: planning.Plan, state: State, jobs: dict[str, str]) -> list[str]:
    """The queued or running jobs whose attempts are of runs that have changed since."""
    definitions = {run.name: run.definition for run in plan.runs}
    stale = set()
    for name, attempts in state.data["attempts"].items():
        for attempt in attempts:
            job = attempt.get("job")
            if (
                job
                and jobs.get(job) in ("pending", "running")
                and attempt["definition"] != definitions.get(name)
            ):
                stale.add(job)
    return sorted(stale)


def cmd_submit(args: argparse.Namespace) -> int:
    campaign, directory, state = open_campaign(args)
    scheduler = scheduler_for(campaign.machine)
    if args.no_build:
        binaries = known_binaries(campaign, state.builds)
        missing = [
            build.name for build in campaign.builds if build.name not in binaries
        ]
        if missing:
            raise CommandError(f"not built yet: {', '.join(missing)}; drop --no-build")
    else:
        try:
            binaries = build_all(campaign, directory, state.builds, log=log)
        finally:
            state.save()
    if not args.no_check:
        result = run_checks(
            campaign, directory, state, binaries, submitting=not args.dry_run
        )
        print_findings(result)
        if not result.ok:
            log("the check failed; nothing is submitted")
            return 1
    plan = write_plan(campaign, directory, state)
    for note in plan.notes:
        log(f"note: {note}")
    jobs = job_states(scheduler, state)
    states = run_states(plan, directory, state, jobs)
    runs, notes = open_runs(campaign, plan, state, states)
    for note in notes:
        log(note)
    stale = stale_jobs(plan, state, jobs)
    if stale:
        log(
            f"note: the jobs {', '.join(stale)} still run attempts of an earlier setting or "
            "build, which the results leave out; 'cancel' cancels the jobs of the campaign"
        )
    if not runs:
        log("no run to submit")
        state.save()
        return 0

    stamp = utc_now().replace(":", "").replace("-", "")
    submitted = []
    for job in plan.jobs:
        job_runs = [run for run in runs if run.job == job.name]
        if not job_runs:
            continue
        entries = []
        for run in job_runs:
            attempt = state.add_attempt(run.name, run.definition)
            entries.append((run, attempt))
        script_dir = directory.jobs / ("dry-run" if args.dry_run else "")
        script = script_dir / f"{job.name}-{stamp}.sh"
        write_text_atomically(
            script,
            measurement_script(
                campaign, directory, scheduler, job.name, entries, binaries
            ),
        )
        script.chmod(0o755)
        if args.dry_run:
            log(f"would submit {script} ({len(entries)} runs)")
            continue
        try:
            job_id = scheduler.submit(script, after=[])
        except SchedulerError as error:
            # the attempts of this job are not made
            for run, attempt in entries:
                state.data["attempts"][run.name].remove(attempt)
            state.save()
            raise CommandError(str(error)) from error
        for _, attempt in entries:
            attempt["job"] = job_id
        state.record_job(
            job_id,
            {
                "name": job.name,
                "kind": "runs",
                "script": str(script),
                "runs": [[run.name, attempt["attempt"]] for run, attempt in entries],
                "state": "pending",
            },
        )
        submitted.append(job_id)
        log(f"submitted {job.name} as job {job_id} ({len(entries)} runs)")
        state.save()

    if not args.dry_run:
        state.save()
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    campaign, directory, state = open_campaign(args)
    if not directory.plan.exists():
        raise CommandError("the campaign has no plan yet; run 'plan' or 'submit'")
    scheduler = scheduler_for(campaign.machine)
    plan = write_plan(campaign, directory, state)
    jobs = job_states(scheduler, state)
    states = run_states(plan, directory, state, jobs)
    state.save()
    counts: dict[str, int] = {}
    log(f"{'run':<40} {'attempt':>7}  {'state':<13} job")
    for run in plan.runs:
        attempt, status = states[run.name]
        counts[status] = counts.get(status, 0) + 1
        number = str(attempt["attempt"]) if attempt else "-"
        job = attempt.get("job") if attempt else None
        log(f"{run.name:<40} {number:>7}  {status:<13} {job or '-'}")
    log(", ".join(f"{count} {status}" for status, count in sorted(counts.items())))
    return 0


def cmd_cancel(args: argparse.Namespace) -> int:
    campaign, directory, state = open_campaign(args)
    scheduler = scheduler_for(campaign.machine)
    jobs = job_states(scheduler, state)
    active = [
        job_id
        for job_id, job_state in jobs.items()
        if job_state in ("pending", "running")
    ]
    scheduler.cancel(active)
    state.save()
    log(f"cancelled {len(active)} job(s)" if active else "no job to cancel")
    return 0


def parser() -> argparse.ArgumentParser:
    main = argparse.ArgumentParser(
        prog="scaling.py",
        description="Strong scaling series of SeisSol, from a campaign file. "
        "See docs/strong-scaling.rst.",
    )
    commands = main.add_subparsers(dest="command", required=True)

    def command(name: str, help_text: str, function) -> argparse.ArgumentParser:
        sub = commands.add_parser(name, help=help_text, description=help_text)
        sub.add_argument(
            "campaign", help="the campaign file, or the directory of a campaign"
        )
        sub.set_defaults(function=function)
        return sub

    sub = command("plan", "plan the runs and jobs, and show them", cmd_plan)
    sub.add_argument("--cells", type=int, default=None, help="the cells of the mesh")
    sub = command(
        "check",
        "check the campaign: layouts, tools, and SeisSol --check of every build and variant",
        cmd_check,
    )
    sub.add_argument(
        "--no-seissol", action="store_true", help="do not run SeisSol --check"
    )
    command("build", "build the executables of the campaign", cmd_build)
    sub = command(
        "submit",
        "build, check, and submit the runs that have not completed (again)",
        cmd_submit,
    )
    sub.add_argument(
        "--dry-run", action="store_true", help="write the job scripts only"
    )
    sub.add_argument(
        "--no-build", action="store_true", help="use the executables built before"
    )
    sub.add_argument("--no-check", action="store_true", help="skip the checks")
    command("status", "show the state of every run", cmd_status)
    command(
        "cancel",
        "cancel the jobs of the campaign that are queued or running",
        cmd_cancel,
    )
    return main


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    with contextlib.ExitStack() as stack:
        args.stack = stack
        try:
            return args.function(args)
        except (ConfigError, CommandError, BuildError, SchedulerError) as error:
            print(f"error: {error}", file=sys.stderr)
            return 2
