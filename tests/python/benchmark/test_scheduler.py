# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Slurm: headers, submissions with dependencies, and the states of jobs."""

from pathlib import Path

import pytest
from strongscaling import config
from strongscaling.scheduler import JobSpec, SchedulerError, SlurmScheduler, slurm_state


def machine(**slurm):
    data = {"scheduler": "slurm", "node": {"cores": 8}, "slurm": slurm}
    return config.parse_machine(data, Path("machine.json"))


@pytest.mark.parametrize(
    "text,state",
    [
        ("PENDING", "pending"),
        ("CONFIGURING", "pending"),
        ("RUNNING", "running"),
        ("COMPLETING", "running"),
        ("COMPLETED", "ended"),
        ("CANCELLED by 1234", "ended"),
        ("TIMEOUT", "ended"),
        ("NODE_FAIL", "ended"),
        ("", "ended"),
    ],
)
def test_states(text, state):
    assert slurm_state(text) == state


def spec(nodes=4):
    return JobSpec(
        name="campaign-n0004",
        nodes=nodes,
        time_limit=3600,
        ranks_per_node=2,
        cpus_per_rank=4,
        gpus_per_node=0,
        log=Path("/c/jobs/n0004-%j.out"),
        workdir=Path("/c"),
    )


def test_header():
    lines = SlurmScheduler(machine(account="p", qos="q", extra=["--ear=off"])).header(
        spec()
    )
    assert "#SBATCH --time=01:00:00" in lines
    assert "#SBATCH --output=/c/jobs/n0004-%j.out" in lines
    assert "#SBATCH --account=p" in lines
    assert "#SBATCH --qos=q" in lines
    assert "#SBATCH --ear=off" in lines
    assert "#SBATCH --no-requeue" in lines
    assert not any("partition" in line for line in lines)


def test_partitions_by_size():
    scheduler = SlurmScheduler(
        machine(
            partition=[
                {"max_nodes": 16, "partition": "micro"},
                {"max_nodes": 768, "partition": "general"},
            ],
            cpus_per_task=False,
        )
    )
    assert "#SBATCH --partition=micro" in scheduler.header(spec(16))
    assert "#SBATCH --partition=general" in scheduler.header(spec(17))
    assert not any("cpus-per-task" in line for line in scheduler.header(spec(16)))


def test_submissions(fake_slurm, tmp_path):
    scheduler = SlurmScheduler(machine())
    script = tmp_path / "job.sh"
    script.write_text("#!/bin/bash\n")
    first = scheduler.submit(script, after=[])
    second = scheduler.submit(script, after=[first, "77"], kind="afterany")
    assert (first, second) == ("1001", "1002")
    submissions = (fake_slurm / "submissions").read_text().splitlines()
    assert submissions[0] == f"--parsable {script}"
    assert submissions[1] == f"--parsable --dependency=afterany:1001:77 {script}"


def test_refused_submission(tmp_path, monkeypatch):
    folder = tmp_path / "bin"
    folder.mkdir()
    (folder / "sbatch").write_text("#!/bin/bash\necho 'invalid account' >&2\nexit 1\n")
    (folder / "sbatch").chmod(0o755)
    monkeypatch.setenv("PATH", str(folder))
    with pytest.raises(SchedulerError, match="invalid account"):
        SlurmScheduler(machine()).submit(tmp_path / "job.sh", after=[])


def test_queue_and_accounting(fake_slurm):
    (fake_slurm / "queue").write_text("11|PENDING\n12|RUNNING\n99|RUNNING\n")
    (fake_slurm / "accounting").write_text("13|COMPLETED\n14|CANCELLED by 5\n")
    states = SlurmScheduler(machine()).states(["11", "12", "13", "14", "15"])
    # a job that is neither queued nor accounted for has ended
    assert states == {
        "11": "pending",
        "12": "running",
        "13": "ended",
        "14": "ended",
        "15": "ended",
    }
