LUMI
====

Website: https://www.lumi-supercomputer.eu/

Here, we concern ourselves with running SeisSol on the **LUMI-G** partition; that is, on the GPU partition.

The nodes then consist of:
* 1× AMD Epyc 7A53 (Zen 3) CPU, configured with 4 NUMA domains
* 4× AMD Instinct MI250x GPUs, thus 8 GCDs in total

Due to the 8 GCDs, we will launch SeisSol with 8 processes per node. The architecture settings we will need for SeisSol are
``milan`` for the CPU architecture (optimizing for Zen 3), and ``gfx90a`` for the GPU architecture (targeting the MI250X).
As device backend, we use HIP, and for the SYCL implementation, we use AdaptiveCpp.

1. Installing Modules

2. Compiling SeisSol

3. Running Jobs

Attached is a job script which does the pinning for us.
The pinning on the LUMI nodes needs some special attention, since 8 out of the 64 cores are reserved for the OS (cf. https://lumi-supercomputer.github.io/LUMI-training-materials/User-Updates/Update-202308/lumig-lownoise/ ).

.. code-block:: bash
    #!/usr/bin/env bash
    #SBATCH --job-name=y-seissol-omp   # Job name
    #SBATCH --output=y-seissol-omp.o # Name of stdout output file
    #SBATCH --error=y-seissol-omp.e  # Name of stderr error file
    #SBATCH --partition=standard-g  # Partition (queue) name
    #SBATCH --nodes=128               # Total number of nodes 
    #SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 128 total (16x8)
    #SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
    #SBATCH --time=01:00:00       # Run time (d-hh:mm:ss)
    #SBATCH --mail-type=all         # Send email at begin and end of job
    #SBATCH --account=<your-project>  # Project for billing
    #SBATCH --mail-user=<your-mail>
    #SBATCH --exclusive
    #SBATCH --requeue


    cat << EOF > select_gpu
    #!/bin/bash

    export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
    exec \$*
    EOF

    chmod +x ./select_gpu

    CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
    CPU_BIND="${CPU_BIND},7e0000,7e000000"
    CPU_BIND="${CPU_BIND},7e,7e00"
    CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

    export MPICH_GPU_SUPPORT_ENABLED=1
    export HSA_XNACK=1

    export OMP_NUM_THREADS=3
    export OMP_PLACES="cores(3)"
    export OMP_PROC_BIND=close

    export DEVICE_STACK_MEM_SIZE=4
    export SEISSOL_FREE_CPUS_MASK="52-54,60-62,20-22,28-30,4-6,12-14,36-38,44-46"

    proxy=SeisSol_proxy_${build_mode}_sgfx90a_hip_6_elastic
    exe=SeisSol_${build_mode}_sgfx90a_hip_6_elastic
    file=parameters.par

    srun --cpu-bind=${CPU_BIND} ./select_gpu ./${exe} ./${file}
