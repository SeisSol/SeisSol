# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""A machine profile drafted from a job script that runs SeisSol."""

from pathlib import Path

from conftest import dedent
from strongscaling import config
from strongscaling.profile import draft_profile, placeholders

LUMI = dedent(
    """
    #!/usr/bin/env bash
    #SBATCH --job-name=seissol   # Job name
    #SBATCH --nodes=4
    #SBATCH --account=project_465000000
    #SBATCH --time=01:00:00
    #SBATCH --partition=standard-g
    #SBATCH --ntasks-per-node=8
    #SBATCH --gpus-per-node=8
    #SBATCH --exclusive
    #SBATCH --requeue

    CPU_BIND="7e000000000000,7e00000000000000"
    CPU_BIND="${CPU_BIND},7e0000,7e000000"

    export MPICH_GPU_SUPPORT_ENABLED=1
    export OMP_NUM_THREADS=3
    export OMP_PLACES="cores(3)"

    srun --unbuffered --cpu-bind=mask_cpu:${CPU_BIND} seissol-launch SeisSol_Release_hip_6_elastic parameters.par
    """
)

SUPERMUC = dedent(
    """
    #!/bin/bash
    #SBATCH -J tpv5
    #SBATCH -o ./%j.%x.out
    #SBATCH --account=pn12ab
    #SBATCH --partition=general
    #SBATCH --nodes=40
    #SBATCH --ntasks-per-node=1
    #SBATCH --ear=off
    module load slurm_setup
    export OMP_NUM_THREADS=94
    export OMP_PLACES="cores(47)"
    ulimit -Ss 2097152
    unset KMP_AFFINITY
    srun SeisSol_Release_sskx_4_elastic parameters.par
    """
)


def test_lumi(tmp_path):
    script = tmp_path / "lumi.sh"
    script.write_text(LUMI)
    profile = draft_profile(script, "lumi")
    assert profile["layout"]["ranks_per_node"] == 8
    assert profile["layout"]["threads_per_rank"] == 3
    assert profile["layout"]["gpus_per_rank"] == 1
    assert profile["node"]["gpus"] == 8
    assert profile["slurm"] == {
        "exclusive": True,
        "account": "project_465000000",
        "partition": "standard-g",
    }
    assert profile["setup"][0] == 'CPU_BIND="7e000000000000,7e00000000000000"'
    assert profile["env"] == {
        "MPICH_GPU_SUPPORT_ENABLED": "1",
        "OMP_PLACES": "cores(3)",
    }
    assert profile["wrapper"] == "seissol-launch"
    assert profile["launcher"] == (
        "srun --nodes={nodes} --ntasks={ranks} --ntasks-per-node={ranks_per_node} "
        "--unbuffered --cpu-bind=mask_cpu:${CPU_BIND}"
    )
    # and a profile that the tool reads
    config.parse_machine(profile, Path("lumi.json"))


def test_supermuc(tmp_path):
    script = tmp_path / "supermuc.sh"
    script.write_text(SUPERMUC)
    profile = draft_profile(script)
    assert profile["layout"]["threads_per_rank"] == 94
    assert profile["slurm"]["extra"] == ["--ear=off"]
    assert profile["setup"] == [
        "module load slurm_setup",
        "ulimit -Ss 2097152",
        "unset KMP_AFFINITY",
    ]
    assert "wrapper" not in profile
    assert profile["launcher"].startswith("srun --nodes={nodes} --ntasks={ranks}")


def test_counts_become_placeholders():
    assert (
        placeholders("mpirun -np 64 --bind-to core")
        == "mpirun -n {ranks} --bind-to core"
    )
    assert placeholders("srun -N 2 -n 16 -c 8") == (
        "srun --ntasks-per-node={ranks_per_node} --nodes={nodes} --ntasks={ranks} "
        "--cpus-per-task={cpus_per_rank}"
    )
    assert placeholders("mpiexec ./x") == "mpiexec -n {ranks} ./x"
