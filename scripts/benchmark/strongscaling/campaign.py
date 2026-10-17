# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The campaign directory, and what it records of the campaign.

    <campaign>/
      campaign.json      the campaign file as read, for the commands run in jobs
      plan.json          the runs and jobs
      state.json         the builds, checks, jobs and attempts so far
      build/ install/    the builds of compile_variants.py
      checks/            the checks of the runs, per build and variant
      jobs/              the job scripts and their logs
      runs/<run>/attempt-<k>/
                         an attempt at a run: its log, its check, its output
      results/           the tables and plots of the campaign

Every command can be repeated: what is recorded is not done again, and a run
is attempted again only if its last attempt did not complete.
"""

from __future__ import annotations

import contextlib
import json
import os
import tempfile
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import fcntl
except ImportError:  # not on Windows
    fcntl = None

STATE_SCHEMA = "seissol-scaling-state/1"

#: the states of an attempt that are final
COMPLETED = "completed"
FAILED_STATES = ("failed", "check-failed", "incomplete", "killed", "not-run")
ACTIVE_STATES = ("queued", "waiting", "running")


@contextlib.contextmanager
def locked(root: Path):
    """Holds the lock of the campaign in @p root, so that no other command changes it meanwhile."""
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".lock").open("w") as stream:
        if fcntl is not None:
            try:
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                print(
                    f"waiting for another command on the campaign in {root}", flush=True
                )
                fcntl.flock(stream, fcntl.LOCK_EX)
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(stream, fcntl.LOCK_UN)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_text_atomically(path: Path, text: str) -> None:
    """Replaces @p path by a file holding @p text in a single step."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            stream.write(text)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def write_json(path: Path, data: Any) -> None:
    write_text_atomically(path, json.dumps(data, indent=2, sort_keys=False) + "\n")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


@dataclass
class Layout:
    root: Path

    @property
    def config_copy(self) -> Path:
        return self.root / "campaign.json"

    @property
    def plan(self) -> Path:
        return self.root / "plan.json"

    @property
    def state(self) -> Path:
        return self.root / "state.json"

    @property
    def builds(self) -> Path:
        return self.root / "build"

    @property
    def installs(self) -> Path:
        return self.root / "install"

    @property
    def checks(self) -> Path:
        return self.root / "checks"

    @property
    def jobs(self) -> Path:
        return self.root / "jobs"

    @property
    def runs(self) -> Path:
        return self.root / "runs"

    @property
    def results(self) -> Path:
        return self.root / "results"

    def attempt(self, run: str, attempt: int) -> Path:
        return self.runs / run / f"attempt-{attempt}"


class State:
    """What the campaign has done so far: builds, checks, jobs and attempts at runs."""

    def __init__(self, path: Path):
        self.path = path
        if path.exists():
            self.data = read_json(path)
            if self.data.get("schema") != STATE_SCHEMA:
                raise ValueError(f"{path} is no state of a campaign of this tool")
        else:
            self.data = {
                "schema": STATE_SCHEMA,
                "builds": {},
                "checks": {},
                "jobs": {},
                "attempts": {},
            }

    def save(self) -> None:
        self.data["updated"] = utc_now()
        write_json(self.path, self.data)

    @property
    def builds(self) -> dict[str, Any]:
        return self.data["builds"]

    @property
    def checks(self) -> dict[str, Any]:
        return self.data["checks"]

    @property
    def jobs(self) -> dict[str, Any]:
        return self.data["jobs"]

    def attempts(self, run: str) -> list[dict[str, Any]]:
        return self.data["attempts"].get(run, [])

    def add_attempt(self, run: str, definition: str) -> dict[str, Any]:
        attempts = self.data["attempts"].setdefault(run, [])
        attempt = {
            "attempt": len(attempts) + 1,
            "runId": str(uuid.uuid4()),
            "definition": definition,
            "job": None,
            "created": utc_now(),
        }
        attempts.append(attempt)
        return attempt

    def record_job(self, job_id: str, info: dict[str, Any]) -> None:
        self.jobs[job_id] = dict(info, submitted=utc_now())


def attempt_status(directory: Path, job_state: str | None) -> str:
    """The state of the attempt in @p directory, of a job in @p job_state.

    job_state is pending, running or ended, or None if it is not known: then an
    attempt that has not finished is unknown.
    """
    exitcode = directory / "exitcode"
    if exitcode.exists():
        if (directory / "check-failed").exists():
            return "check-failed"
        try:
            code = int(exitcode.read_text(encoding="utf-8").strip() or "1")
        except ValueError:
            code = 1
        if code != 0:
            return "failed"
        summary = directory / "output" / "data-summary.json"
        try:
            status = read_json(summary).get("status")
        except (OSError, ValueError):
            status = None
        return COMPLETED if status == "completed" else "incomplete"
    if (directory / "started").exists():
        if job_state == "running":
            return "running"
        return "killed" if job_state == "ended" else "unknown"
    if job_state == "pending":
        return "queued"
    if job_state == "running":
        return "waiting"
    if job_state == "ended":
        return "not-run"
    return "unknown" if directory.exists() else "planned"
