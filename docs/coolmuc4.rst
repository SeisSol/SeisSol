..
  SPDX-FileCopyrightText: 2019 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_coolmuc4:


CoolMUC4
===========

CoolMUC4 offers 106 Intel(R) Xeon(R) Platinum 8480+ compute nodes with 112 cores per node or 12 Intel(R) Xeon(R) Platinum 8380 CPU nodes with 80 cores per node.
As the precompiled modules have mismatch in compiler versions quite often, we build all the dependencies using spack using a particular compiler from scratch in this installation workflow.

Note: This is one way to install SeisSol on CoolMUC4 with spack. One could also install all dependencies in a shared location manually.

Installing Git
---------------------------------------
Even though Git is available on CoolMUC4, it is not available on the compute node if you want to load the spack environment there.
Therefore, we recommend installing Git in a separate location, for example in a folder called ``$HOME/seissol_env/local``.
.. code-block:: bash

    mkdir -p $HOME/seissol_env/local

    export SEISSOL_BASE=$HOME/seissol_env
    export SEISSOL_PREFIX=$SEISSOL_BASE/local
    export PATH=$SEISSOL_PREFIX/bin:$PATH

    wget https://mirrors.edge.kernel.org/pub/software/scm/git/git-2.43.0.tar.gz
    tar -xvzf git-2.43.0.tar.gz
    cd git-2.43.0
    make prefix=$SEISSOL_PREFIX all
    make prefix=$SEISSOL_PREFIX install
    cd

Setting up the Spack environment
---------------------------------------
We recommend installing SeisSol and its dependencies using Spack.
First, clone Spack from its Git repository:
.. code-block:: bash

    git clone https://github.com/spack/spack.git
    . ~/spack/share/spack/setup-env.sh
    module load intel/2024.1.0
    module load gcc/13.2.0
    spack compiler find

Next, create a spack environment for SeisSol with ``spack env create seissol``. Then activate the environment with ``spacktivate -p seissol``.
Next, add the spack requirements for SeisSol by running ``spack config edit`` and adding the following packages to specs

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
Once installed, change the environment variable ``CC`` and load intel-oneapi-mpi using:
.. code-block:: bash
  export CC=icx
  export CXX=icpx
  spack load intel-oneapi-mpi

Compiling SeisSol
---------------------------------------
Once all dependencies are installed, you can compile SeisSol with the following steps:
.. code-block:: bash
    git clone https://github.com/SeisSol/SeisSol.git
    cd SeisSol
    git submodule update --init --recursive
    mkdir -p build
    cd build
    cmake .. -DHOST_ARCH=skx -DCMAKE_BUILD_TYPE=Release -DORDER=4 -DEQUATIONS=elastic
    make -j

This will create the SeisSol executable ``SeisSol_Release_skx_4_elastic`` in the build directory.

Running SeisSol
---------------

It is recommended to add the following lines to your ``~/.bashrc`` to load the spack environment and the required modules automatically 
once the entire installation is done. This would load the required modules when you log in:

.. code-block:: bash

  export SEISSOL_BASE=$HOME/seissol_env
  export SEISSOL_PREFIX=$SEISSOL_BASE/local
  export PATH=$SEISSOL_PREFIX/bin:$PATH

  . ~/spack/share/spack/setup-env.sh
  module load intel/2024.1.0
  module load gcc/13.2.0

  spacktivate -p seissol
  export CC=icx
  export CXX=icpx
  spack load intel-oneapi-mpi

This is an example job submission script for SeisSol on CoolMUC4. For your applications, change ``#SBATCH --nodes=``
to the number of nodes you want to run on. A rule of thumb for optimal performance is to distribute your jobs to 1 node per 100k elements. This rule of thumb does not account for potentially shorter queue times.

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

  #Run the program:
  export I_MPI_PIN=0 # disable automatic pinning by I_MPI to allow OMP_PLACES to work correctly
  export OMP_NUM_THREADS=110
  export OMP_PLACES="cores(55)"
  export OMP_PROC_BIND="close"
  module load slurm_setup
  unset KMP_AFFINITY

  echo 'num_nodes:' $SLURM_JOB_NUM_NODES 'ntasks:' $SLURM_NTASKS

  SEISSOL_EXECUTABLE=~/SeisSol/build/SeisSol_Release_skx_4_elastic

  mpirun -n $SLURM_NTASKS $SEISSOL_EXECUTABLE parameters.par
