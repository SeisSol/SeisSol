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

    assert cli.main(["report", str(path)]) == 0
    results = path.parent / "fake" / "results"
    runs = rows(results / "runs.csv")
    assert len(runs) == 12
    assert {run["status"] for run in runs} == {"completed"}
    assert {run["cellsPerRank"] for run in runs if run["nodes"] == "4"} == {"250.0"}
    scaling = rows(results / "scaling.csv")
    four = [
        entry for entry in scaling if entry["nodes"] == "4" and entry["variant"] == "a"
    ][0]
    # the fake takes 100 / ranks + 2 seconds
    assert float(four["time"]) == 27.0
    assert abs(float(four["speedup"]) - 102.0 / 27.0) < 1e-9
    assert len(rows(results / "ranks.csv")) == 2 * 2 * (1 + 2 + 4)
    assert len(rows(results / "epochs.csv")) == 12 * 4
    report = (results / "report.md").read_text()
    assert "## fake/a" in report
    summary = json.loads((results / "scaling.json").read_text())
    assert summary["checks"] == []


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
        job={"time_limit": 30, "report_job": True},
    )
    assert cli.main(["submit", str(path)]) == 0
    out = capsys.readouterr().out
    assert "submitted n0004 as job 1003" in out
    assert "submitted the report as job 1004" in out
    submissions = (fake_slurm / "submissions").read_text().splitlines()
    assert submissions[-1].startswith(
        "--parsable --dependency=afterany:1001:1002:1003 "
    )

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
    assert "submitted n0001 as job 1005" in capsys.readouterr().out
    (fake_slurm / "queue").write_text("1005|PENDING\n1006|PENDING\n")
    assert cli.main(["cancel", str(path)]) == 0
    assert "cancelled 2 job(s)" in capsys.readouterr().out
    assert (fake_slurm / "cancelled").read_text().split() == ["1005", "1006"]


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


def test_a_local_report_job(campaign_file):
    path = campaign_file(
        scaling={"nodes": [1, 2], "repetitions": 1},
        variants=[{"name": "a"}],
        job={"time_limit": 10, "report_job": True},
    )
    assert cli.main(["submit", str(path)]) == 0
    assert (path.parent / "fake" / "results" / "report.md").exists()
