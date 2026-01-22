..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_coolmuc4:

CoolMUC4
========

CoolMUC4 offers 106 Intel(R) Xeon(R) Platinum 8480+ compute nodes with 112 cores per node with two-way hyperthreading and 12 Intel(R) Xeon(R) Platinum 8380 CPU nodes with 80 cores per node. The official documentation can be found at: https://doku.lrz.de/coolmuc-4-1082337877.html

As the provided precompiled modules have mismatch in compiler versions quite often, we build all the dependencies using spack using a particular compiler from scratch in this installation workflow.

Note: This is one way to install SeisSol on CoolMUC4 with spack. One could also install all dependencies in a private location manually.

Installing Git
--------------
Even though Git is available on CoolMUC4, it is not available on the compute node if you want to load the spack environment there.
Therefore, we recommend installing Git in a separate location, for example in a folder called ``$HOME/seissol-env/local``.

.. code-block:: bash

    export SEISSOL_BASE=$HOME/seissol-env
    export SEISSOL_PREFIX=$SEISSOL_BASE/local
    export PATH=$SEISSOL_PREFIX/bin:$PATH
    mkdir -p $SEISSOL_BASE
    mkdir -p $SEISSOL_PREFIX

    wget https://mirrors.edge.kernel.org/pub/software/scm/git/git-2.43.0.tar.gz
    tar -xvzf git-2.43.0.tar.gz
    cd git-2.43.0
    make prefix=$SEISSOL_PREFIX all
    make prefix=$SEISSOL_PREFIX install
    cd

Setting up the Spack environment
--------------------------------
We recommend installing dependencies of SeisSol using Spack.
First, clone Spack from its Git repository:

.. code-block:: bash

    git clone https://github.com/spack/spack.git
    . <where spack was installed>/spack/share/spack/setup-env.sh
    module load intel/2024.1.0
    module load gcc/13.2.0
    spack compiler find

Next, create a spack environment for SeisSol with ``spack env create seissol``. Then activate the environment with ``spacktivate -p seissol``.
Next, add the spack requirements for SeisSol by running ``spack config edit`` and adding the following packages to specs:

.. code-block:: yaml

    - cmake
    - ninja
    - easi%oneapi ^intel-oneapi-mpi ^asagi~fortran
    - hdf5 ^intel-oneapi-mpi
    - libxsmm+generator%oneapi
    - netcdf-c ^intel-oneapi-mpi
    - eigen
    - metis +int64
    - parmetis +int64 ^intel-oneapi-mpi
    - python@3.12%gcc@13.2.0
    - py-numpy%gcc@13.2.0
    - py-scipy%gcc@13.2.0
    - py-setuptools
    - re2c
    - bison
    - py-pspamm

Run ``spack concretize -f`` and ``spack install`` to install all dependencies. This might take a while, so be patient.
Once installed, change the environment variable ``CC``, ``CXX`` and ``FC`` and load intel-oneapi-mpi using:

.. code-block:: bash

    export CC=icx
    export CXX=icpx
    export FC=ifx
    spack load intel-oneapi-mpi

Compiling SeisSol
-----------------
Once all dependencies are installed, you can compile SeisSol with the following steps:

.. code-block:: bash

    git clone https://github.com/SeisSol/SeisSol.git
    cd SeisSol
    git submodule update --init --recursive
    mkdir -p build
    cd build
    cmake .. -DHOST_ARCH=skx -DCMAKE_BUILD_TYPE=Release -DORDER=4 -DEQUATIONS=elastic
    make -j

This will create the SeisSol executable according to CMake flag ``NEW_BINARY_NAMING`` you chose in the build directory. The default ``OFF`` option will create the executable with the name ``SeisSol_Release_skx_4_elastic``. If you choose ``ON``, the executable will have the name ``seissol-cpu-elastic-p4-f64``.

Running SeisSol
---------------
When working with SeisSol, we recommend to create an extra file in ``SEISSOL_BASE`` with the following content and source it in ``~/.bashrc`` to load the spack environment and the required modules automatically once the entire installation is done. This would load the required modules when you log in:

.. code-block:: bash

    export SEISSOL_BASE=$HOME/seissol-env
    export SEISSOL_PREFIX=$SEISSOL_BASE/local
    export PATH=$SEISSOL_PREFIX/bin:$PATH

    . <where spack is installed>/spack/share/spack/setup-env.sh
    module load intel/2024.1.0
    module load gcc/13.2.0

    spacktivate -p seissol
    export CC=icx
    export CXX=icpx
    spack load intel-oneapi-mpi

For quick testing, one could request for a node via ``salloc -M inter -p cm4_inter --nodes=2 --ntasks-per-node=2 --cpus-per-task=112 -t 00:30:00`` and run SeisSol directly on the compute node using

.. code-block:: bash

    LOGICAL_CORES=$(expr $SLURM_CPUS_PER_TASK - 1)
    PHYSICAL_CORES=$((NUM_COMPUTE_CORES)/2)
    export I_MPI_PIN=0 # disable automatic pinning by I_MPI to allow OMP_PLACES to work correctly
    export OMP_NUM_THREADS=$LOGICAL_CORES
    export OMP_PLACES="cores($PHYSICAL_CORES)"
    export OMP_PROC_BIND="close"
    unset KMP_AFFINITY

    SEISSOL_EXECUTABLE=~/SeisSol/build/SeisSol_Release_skx_4_elastic

    mpirun -n $SLURM_NTASKS $SEISSOL_EXECUTABLE parameters.par

For bigger runs, this is an example job submission script for SeisSol on CoolMUC4. For your applications, change ``#SBATCH --nodes=`` to the number of nodes you want to run on. A rule of thumb for optimal performance is to distribute your jobs to 1 node per 100k elements.

.. code-block:: bash

    #!/bin/bash
    # Job Name and Files (also --job-name)
    #SBATCH -J <job name>

    #Output and error (also --output, --error):
    #SBATCH -o ./%j.%x.out
    #SBATCH -e ./%j.%x.err

    #Initial working directory:
    #SBATCH --chdir=<work directory>

    #Notification and type
    #SBATCH --mail-type=END
    #SBATCH --mail-user=<your email address>

    #Number of nodes and MPI tasks per node:
    #SBATCH --clusters=cm4
    #SBATCH --partition=cm4_std
    #SBATCH --qos=cm4_std
    #SBATCH --nodes=2
    #SBATCH --ntasks-per-node=2
    #SBATCH --cpus-per-task=112
    #SBATCH --time=03:00:00

    source ~/.bashrc
    module list

    #Setup pinning:
    LOGICAL_CORES=$(expr $SLURM_CPUS_PER_TASK - 1)
    PHYSICAL_CORES=$((NUM_COMPUTE_CORES)/2)
    export I_MPI_PIN=0 # disable automatic pinning by I_MPI to allow OMP_PLACES to work correctly
    export OMP_NUM_THREADS=$LOGICAL_CORES
    export OMP_PLACES="cores($PHYSICAL_CORES)"
    export OMP_PROC_BIND="close"
    unset KMP_AFFINITY

    echo 'num_nodes:' $SLURM_JOB_NUM_NODES 'ntasks:' $SLURM_NTASKS

    #Run SeisSol:
    SEISSOL_EXECUTABLE=~/SeisSol/build/SeisSol_Release_skx_4_elastic

    mpirun -n $SLURM_NTASKS $SEISSOL_EXECUTABLE parameters.par
