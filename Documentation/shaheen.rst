..
  SPDX-FileCopyrightText: 2022-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_shaheen:


Shaheen
=======


Add these lines to ``~/.bashrc``:

::

    ##### module load for SeisSol
    module unload PrgEnv-cray
    module load PrgEnv-gnu
    module load python
    module load cmake
    module unload intel
    # to use preinstall seissol-env module
    module use /project/k1587/ulrich/spack/share/spack/modules/cray-cnl7-ivybridge/
    module load seissol-env-develop-gcc-11.2.0-rckyrcj
    export CC=cc
    export CXX=CC
    export FC=ftn

This will load a preinstalled seissol-env module.

Alternatively (and for reference), to compile seissol-env on shaheen, follow the procedure below:

.. code-block:: bash

    git clone --depth 1 --branch v0.18.1 https://github.com/spack/spack.git
    cd spack
    echo "export SPACK_ROOT=$PWD" >> $HOME/.bashrc
    echo "export PATH=\$SPACK_ROOT/bin:\$PATH" >> $HOME/.bashrc
    # clone seissol-spack-aid and add the repository
    git clone --branch supermuc_NG https://github.com/SeisSol/seissol-spack-aid.git
    cd seissol-spack-aid
    spack repo add ./spack
    spack compiler find


Then update ``~/.spack/packages.yaml`` as follow:

.. code-block:: yaml

    packages:
      python:
        externals:
        - spec: python@3.10.1
          buildable: False
          modules:
           - python/3.10.1-cdl
      cmake:
        buildable: false
        externals:
        - spec: cmake@3.22.1
          modules:
           - cmake/3.22.1
      mpich:
        buildable: false
        externals:
        - spec: mpich@7.7.18
          modules:
          - cray-mpich/7.7.18
        - spec: mpich@7.7.20
          modules:
          - cray-mpich/7.7.20
      all:
        providers:
          mpi: [mpich]


Finally, install seissol-env with

.. code-block:: bash

    spack install -j 8 seissol-env %gcc@11.2.0 ^mpich

and create a module with:

.. code-block:: bash

    spack module tcl refresh seissol-env@develop%%gcc@11.2.0

To access the module at start up, add to your ``~/.bashrc``:

.. code-block:: bash

    module use $SPACK_ROOT/share/spack/modules/cray-cnl7-ivybridge/

Finally, install SeisSol with cmake, as usual, with ``-DHOST_ARCH=hsw``.

Here is an example job submission script for SeisSol on Shaheen (to be launched from the ``/scratch/`` folder):

::

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

  # Wall clock limit:
  #SBATCH --time=00:30:00
  #SBATCH --no-requeue

  #Setup of execution environment
  #SBATCH --export=ALL
  #SBATCH --account=<project id>
  #SBATCH --partition=debug

  #Number of nodes and MPI tasks per node:
  #SBATCH --nodes=2
  #SBATCH --ntasks-per-node=1

  export MP_SINGLE_THREAD=no
  unset KMP_AFFINITY
  export OMP_NUM_THREADS=31
  # you could also consider OMP_NUM_THREADS=62 for high order
  export OMP_PLACES="cores(31)"
  #Prevents errors such as experience in Issue #691
  export I_MPI_SHM_HEAP_VSIZE=8192

  export XDMFWRITER_ALIGNMENT=8388608
  export XDMFWRITER_BLOCK_SIZE=8388608
  export SC_CHECKPOINT_ALIGNMENT=8388608

  export SEISSOL_CHECKPOINT_ALIGNMENT=8388608
  export SEISSOL_CHECKPOINT_DIRECT=1
  export ASYNC_MODE=THREAD
  export ASYNC_BUFFER_ALIGNMENT=8388608

  export MPICH_MAX_THREAD_SAFETY=multiple
  # update to output folder
  lfs setstripe -c 32 output
  srun path_2_SeisSol_Release_dhsw_4_elastic parameters.par
