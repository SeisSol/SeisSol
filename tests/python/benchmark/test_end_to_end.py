# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""A whole campaign of the fake SeisSol on the local machine, and on a fake Slurm."""

import csv
import json

from strongscaling import cli


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_a_local_campaign(campaign_file, capsys):
    path = campaign_file()
    assert cli.main(["plan", str(path)]) == 0
    assert "12 runs in 3 jobs" in capsys.readouterr().out
    assert cli.main(["check", str(path)]) == 0
    out = capsys.readouterr().out
    assert "the check passed" in out
    assert "1000 cells" in out

    assert cli.main(["submit", str(path)]) == 0
    assert cli.main(["status", str(path)]) == 0
    assert "12 completed" in capsys.readouterr().out
    # submitting again submits nothing
    assert cli.main(["submit", str(path), "--no-check"]) == 0
    assert "no run to submit" in capsys.readouterr().out


def test_failed_runs_are_attempted_again(campaign_file, tmp_path, monkeypatch, capsys):
    path = campaign_file(
        scaling={"nodes": [1], "repetitions": 1}, variants=[{"name": "a"}]
    )
    monkeypatch.setenv("FAKE_FAIL_ONCE", str(tmp_path / "failed-once"))
    assert cli.main(["submit", str(path)]) == 0
    assert cli.main(["status", str(path)]) == 0
    assert "1 failed" in capsys.readouterr().out
    assert cli.main(["submit", str(path)]) == 0
    assert cli.main(["status", str(path)]) == 0
    out = capsys.readouterr().out
    assert "1 completed" in out
    assert "n0001-fake-a-r1" in out and "      2  completed" in out


def test_a_failing_check_submits_nothing(campaign_file, monkeypatch, capsys):
    path = campaign_file()
    monkeypatch.setenv("FAKE_CHECK_FAILS", "1")
    assert cli.main(["submit", str(path)]) == 1
    out = capsys.readouterr().out
    assert "materialfilename" in out
    assert "nothing is submitted" in out
    assert not (path.parent / "fake" / "runs").exists()


def test_a_slurm_campaign(campaign_file, fake_slurm, capsys):
    path = campaign_file(
        machine={
            "scheduler": "slurm",
            "node": {"cores": 4},
            "layout": {"ranks_per_node": 1, "cpus_per_rank": 4},
            "slurm": {"account": "project"},
        },
        job={"time_limit": 30},
    )
    assert cli.main(["submit", str(path)]) == 0
    out = capsys.readouterr().out
    assert "submitted n0004 as job 1003" in out

    # queued: nothing is submitted again
    (fake_slurm / "queue").write_text("1001|PENDING\n1002|PENDING\n1003|RUNNING\n")
    assert cli.main(["submit", str(path), "--no-check"]) == 0
    assert "no run to submit" in capsys.readouterr().out
    assert cli.main(["status", str(path)]) == 0
    out = capsys.readouterr().out
    assert "8 queued, 4 waiting" in out

    # the jobs ended without running: every run is attempted again
    (fake_slurm / "queue").write_text("")
    (fake_slurm / "accounting").write_text(
        "1001|CANCELLED\n1002|TIMEOUT\n1003|FAILED\n"
    )
    assert cli.main(["status", str(path)]) == 0
    assert "12 not-run" in capsys.readouterr().out
    assert cli.main(["submit", str(path), "--no-check"]) == 0
    assert "submitted n0001 as job 1004" in capsys.readouterr().out
    (fake_slurm / "queue").write_text("1004|PENDING\n1005|PENDING\n")
    assert cli.main(["cancel", str(path)]) == 0
    assert "cancelled 2 job(s)" in capsys.readouterr().out
    assert (fake_slurm / "cancelled").read_text().split() == ["1004", "1005"]


def test_a_dry_run(campaign_file, capsys):
    path = campaign_file(
        machine={
            "scheduler": "slurm",
            "node": {"cores": 4},
            "layout": {"ranks_per_node": 1, "cpus_per_rank": 4},
        }
    )
    assert cli.main(["submit", str(path), "--dry-run"]) == 0
    assert "would submit" in capsys.readouterr().out
    scripts = sorted((path.parent / "fake" / "jobs" / "dry-run").glob("*.sh"))
    assert len(scripts) == 3
    state = json.loads((path.parent / "fake" / "state.json").read_text())
    assert state["attempts"] == {}


def test_another_executable_is_another_series(campaign_file, fakes, capsys):
    path = campaign_file(
        scaling={"nodes": [1, 2], "repetitions": 1}, variants=[{"name": "a"}]
    )
    assert cli.main(["submit", str(path)]) == 0
    assert cli.main(["submit", str(path)]) == 0
    assert "no run to submit" in capsys.readouterr().out
    # rebuilt: the runs of the earlier executable do not count any more
    binary = fakes["seissol"]
    binary.write_text(binary.read_text() + "\n# rebuilt\n")
    assert cli.main(["status", str(path)]) == 0
    assert "2 planned" in capsys.readouterr().out
    assert cli.main(["submit", str(path)]) == 0
    assert cli.main(["status", str(path)]) == 0
    assert "2 completed" in capsys.readouterr().out


def test_an_unreadable_queue_submits_nothing(campaign_file, fake_slurm, capsys):
    path = campaign_file(
        machine={
            "scheduler": "slurm",
            "node": {"cores": 4},
            "layout": {"ranks_per_node": 1, "cpus_per_rank": 4},
        },
        scaling={"nodes": [1], "repetitions": 1},
        variants=[{"name": "a"}],
    )
    assert cli.main(["submit", str(path), "--no-check"]) == 0
    (fake_slurm / "squeue").write_text(
        "#!/usr/bin/env bash\necho 'slurm_load_jobs error' >&2\nexit 1\n"
    )
    assert cli.main(["submit", str(path), "--no-check"]) == 0
    out = capsys.readouterr().out
    assert "squeue failed" in out
    assert "no run to submit" in out
    assert (fake_slurm / "count").read_text().strip() == "1"


def test_an_executable_for_other_nodes(campaign_file, monkeypatch, capsys):
    path = campaign_file()
    monkeypatch.setenv("FAKE_SIGNAL", "1")
    assert cli.main(["check", str(path)]) == 0
    out = capsys.readouterr().out
    assert "does not run on this node (SIGILL)" in out
    assert "the check passed" in out
