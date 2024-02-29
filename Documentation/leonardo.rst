Leonardo
========

Website: https://leonardo-supercomputer.cineca.eu/

We will focus here on the Booster module of Leonardo, i.e. we use GPUs. That is, our nodes consist of:
* 1× Intel SPR???? (Sapphire Rapids) CPU
* 4× Nvidia A100 SXM6 GPUs

Thus, we will run SeisSol with 4 ranks per node. As architectures, we compile for the host/CPU architecture ``skx``, and use ``sm_80`` as architecture for the GPUs, together
with CUDA as device backend. For the SYCL parts, we use AdaptiveCpp (formerly known as hipSYCL or Open SYCL).

1. Installing Modules

TODO

2. Compiling SeisSol

TODO

3. Running Jobs

Attached is a job script which does the pinning for us.

.. code-block:: bash
    #!/usr/bin/env bash
    #SBATCH --account=<PROJECT_NAME>
    #SBATCH --qos=normal
    #SBATCH --partition=boost_usr_prod
    #SBATCH --time 01:00:00
    #SBATCH --nodes=8
    #SBATCH --ntasks-per-node=4
    #SBATCH --cpus-per-task=8
    #SBATCH --gres=gpu:4
    #SBATCH --mem=200G
    #SBATCH --job-name=seissol-gpu
    #SBATCH --exclusive
    #SBATCH --output=y-seissol.o
    #SBATCH --error=y-seissol.e
    #SBATCH --export=ALL

    export OMP_NUM_THREADS=4
    export OMP_PLACES="cores(4)"
    export OMP_BIND="spread"

    export OMP_PROC_BIND=close
    export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
    export SLURM_CPU_BIND_TYPE="cores"

    export DEVICE_STACK_MEM_SIZE=2
    export SEISSOL_FREE_CPUS_MASK="16-19,20-23,24-27,28-31"

    cat << EOF > select_gpu
    #!/bin/bash

    export CUDA_VISIBLE_DEVICES=\$SLURM_LOCALID
    "\$@"
    EOF

    chmod +x ./select_gpu

    proxy=$(basename $(realpath ./SeisSol_p*))
    exe=$(basename $(realpath ./SeisSol_[^p]*))
    file=parameters.par

    ulimit -c unlimited

    CPU_BIND="mask_cpu"
    CPU_BIND="${CPU_BIND}:000f000f"
    CPU_BIND="${CPU_BIND},00f000f0"
    CPU_BIND="${CPU_BIND},0f000f00"
    CPU_BIND="${CPU_BIND},f000f000"

    srun --cpu-bind=${CPU_BIND} ./select_gpu ./${exe} ./${file}

