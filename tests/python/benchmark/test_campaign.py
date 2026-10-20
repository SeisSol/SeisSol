# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The attempts at runs: their states, from what they left, and their records."""

import json

from strongscaling.campaign import State, attempt_status


def test_attempt_states(tmp_path):
    folder = tmp_path / "attempt-1"
    assert attempt_status(folder, None) == "planned"
    assert attempt_status(folder, "pending") == "queued"
    assert attempt_status(folder, "running") == "waiting"
    assert attempt_status(folder, "ended") == "not-run"
    folder.mkdir()
    (folder / "started").write_text("now\n")
    assert attempt_status(folder, "running") == "running"
    assert attempt_status(folder, "ended") == "killed"
    assert attempt_status(folder, None) == "unknown"
    (folder / "exitcode").write_text("0\n")
    assert attempt_status(folder, "ended") == "incomplete"
    (folder / "output").mkdir()
    (folder / "output" / "data-summary.json").write_text(
        json.dumps({"status": "completed"})
    )
    assert attempt_status(folder, "ended") == "completed"
    (folder / "exitcode").write_text("134\n")
    assert attempt_status(folder, "ended") == "failed"
    (folder / "check-failed").write_text("1\n")
    assert attempt_status(folder, "ended") == "check-failed"


def test_attempts_are_recorded(tmp_path):
    state = State(tmp_path / "state.json")
    first = state.add_attempt("run", "abc")
    second = state.add_attempt("run", "abc")
    assert (first["attempt"], second["attempt"]) == (1, 2)
    assert first["runId"] != second["runId"]
    state.save()
    assert State(tmp_path / "state.json").attempts("run")[1]["runId"] == second["runId"]
