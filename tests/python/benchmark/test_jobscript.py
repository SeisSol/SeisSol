# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Job scripts start every run with its layout, environment and overrides."""

import subprocess

from strongscaling import config, jobscript, plan
from strongscaling.campaign import Layout, State
from strongscaling.scheduler import LocalScheduler, SlurmScheduler


def test_launcher_placeholders(campaign_file):
    campaign = config.load_campaign(campaign_file())
    run = plan.make_plan(campaign).runs[-1]
    filled = jobscript.fill_launcher(
        "srun -N {nodes} -n {ranks} --ntasks-per-node={ranks_per_node} "
        "-c {cpus_per_rank} --cpu-bind=mask_cpu:${CPU_BIND} {unknown}",
        run,
    )
    assert filled == (
        "srun -N 4 -n 4 --ntasks-per-node=1 -c 4 --cpu-bind=mask_cpu:${CPU_BIND} {unknown}"
    )


def test_overrides(campaign_file):
    campaign = config.load_campaign(
        campaign_file(
            case={
                "parameter_file": "case/parameters.par",
                "set": {"abortcriteria.endtime": 1.0, "output.wavefieldoutput": 0},
                "epoch_interval": 0.25,
            },
            variants=[{"name": "a", "set": {"abortcriteria.endtime": 2.0}}],
        )
    )
    overrides = jobscript.run_overrides(campaign, campaign.variants[0], "/out/data")
    assert overrides == {
        "abortcriteria.endtime": "2.0",
        "output.wavefieldoutput": "0",
        "output.epochinterval": "0.25",
        "output.outputfile": "/out/data",
    }


def test_time_limits():
    assert (
        jobscript.with_time_limit("srun -n 4", "slurm", 600)
        == "srun --time=00:10:00 -n 4"
    )
    assert (
        jobscript.with_time_limit("mpirun -np 4", "local", 600)
        == "timeout 600 mpirun -np 4"
    )
    assert jobscript.with_time_limit("srun", "slurm", None) == "srun"


def entries(campaign, directory, job_name):
    the_plan = plan.make_plan(campaign)
    state = State(directory.state)
    runs = [run for run in the_plan.runs if run.job == job_name]
    return [(run, state.add_attempt(run.name, run.definition)) for run in runs]


def test_variants_of_other_layouts_share_a_job(campaign_file):
    campaign = config.load_campaign(
        campaign_file(
            machine={
                "scheduler": "slurm",
                "node": {"cores": 8},
                "layout": {"ranks_per_node": 1, "cpus_per_rank": 8},
            },
            variants=[{"name": "one"}, {"name": "four", "ranks_per_node": 4}],
        )
    )
    directory = Layout(campaign.directory)
    text = jobscript.measurement_script(
        campaign,
        directory,
        SlurmScheduler(campaign.machine),
        "n0001",
        entries(campaign, directory, "n0001"),
        {"fake": campaign.builds[0].binary},
    )
    # the layout of the run of four ranks of two CPUs, which fits a node of eight
    assert "#SBATCH --ntasks-per-node=4" in text
    assert "#SBATCH --cpus-per-task=2" in text
    assert "--ntasks-per-node=1 --cpus-per-task=8" in text
    assert "--ntasks-per-node=4 --cpus-per-task=2" in text


def test_a_slurm_job(campaign_file, tmp_path):
    campaign = config.load_campaign(
        campaign_file(
            machine={
                "scheduler": "slurm",
                "node": {"cores": 4, "gpus": 2},
                "layout": {"ranks_per_node": 2, "cpus_per_rank": 2, "gpus_per_rank": 1},
                "slurm": {
                    "account": "project",
                    "partition": "big",
                    "extra": ["--ear=off"],
                },
                "setup": ["module load seissol"],
                "env": {"OMP_PLACES": "cores(1)"},
                "wrapper": "seissol-launch",
            },
            variants=[
                {"name": "a", "env": {"SEISSOL_COMMTHREAD": "0", "QUOTED": "a b"}}
            ],
        )
    )
    directory = Layout(campaign.directory)
    text = jobscript.measurement_script(
        campaign,
        directory,
        SlurmScheduler(campaign.machine),
        "n0002",
        entries(campaign, directory, "n0002"),
        {"fake": campaign.builds[0].binary},
    )
    assert "#SBATCH --nodes=2" in text
    assert "#SBATCH --ntasks-per-node=2" in text
    assert "#SBATCH --cpus-per-task=2" in text
    assert "#SBATCH --gpus-per-node=2" in text
    assert "#SBATCH --account=project" in text
    assert "#SBATCH --partition=big" in text
    assert "#SBATCH --ear=off" in text
    assert "#SBATCH --exclusive" in text
    assert "module load seissol" in text
    assert "export OMP_PLACES='cores(1)'" in text
    assert "export QUOTED='a b'" in text
    assert "export OMP_NUM_THREADS=1" in text
    assert (
        'srun --nodes=2 --ntasks=4 --ntasks-per-node=2 --cpus-per-task=2 seissol-launch "$seissol" "${arguments[@]}"'
        in text
    )
    assert text.count("SEISSOL_RUN_ID=") == 2
    # the script is valid bash
    script = tmp_path / "job.sh"
    script.write_text(text)
    assert subprocess.run(["bash", "-n", str(script)], check=False).returncode == 0


def test_a_local_job_runs(campaign_file, fakes):
    campaign = config.load_campaign(campaign_file())
    directory = Layout(campaign.directory)
    scheduler = LocalScheduler(campaign.machine)
    job = entries(campaign, directory, "n0002")
    text = jobscript.measurement_script(
        campaign, directory, scheduler, "n0002", job, {"fake": fakes["seissol"]}
    )
    script = directory.jobs / "n0002.sh"
    script.parent.mkdir(parents=True)
    script.write_text(text)
    scheduler.submit(script, after=[])
    for run, attempt in job:
        folder = directory.attempt(run.name, attempt["attempt"])
        assert (folder / "exitcode").read_text().strip() == "0"
        assert (folder / "check.json").exists()
        assert (folder / "output" / "data-summary.json").exists()
        assert "on 2 ranks" in (folder / "seissol.log").read_text()


def test_the_report_job(campaign_file):
    campaign = config.load_campaign(campaign_file())
    directory = Layout(campaign.directory)
    text = jobscript.report_script(
        campaign, directory, SlurmScheduler(campaign.machine)
    )
    assert "#SBATCH --nodes=1" in text
    assert "#SBATCH --exclusive" not in text
    assert f"report {directory.root}" in text
