..
  SPDX-FileCopyrightText: 2022-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_frontera:


Frontera
========


Add these lines to ``~/.bashrc``:

::

    module switch python3 python3/3.9.2
    module use /work2/09160/ulrich/frontera/spack/share/spack/modules/linux-centos7-cascadelake
    module load seissol-env
    export CC=mpiicc
    export CXX=mpiicpc

This will load a preinstalled seissol-env module.

Alternatively (and for reference), to compile seissol-env on Frontera, follow the procedure below:

.. code-block:: bash

    git clone --depth 1 --branch v0.21.0 https://github.com/spack/spack.git
    cd spack
    echo "export SPACK_ROOT=$PWD" >> $HOME/.bashrc
    echo "export PATH=\$SPACK_ROOT/bin:\$PATH" >> $HOME/.bashrc
    # clone seissol-spack-aid and add the repository
    git clone https://github.com/SeisSol/seissol-spack-aid.git
    spack repo add seissol-spack-aid/spack
    spack compiler find

Following the workaround proposed in https://github.com/spack/spack/issues/10308, precise the module of the intel compilers in ``~/.spack/linux/compilers.yaml`` by changing ``modules: []`` to ``modules: ['intel/19.1.1']``.

Then, update ``~/.spack/packages.yaml`` as follow:

.. code-block:: yaml

    packages:
      python:
        externals:
        - spec: python@3.9.2
          buildable: false
          prefix: /opt/apps/intel19/python3/3.9.2
          modules:
          - python3/3.9.2
      intel-mpi:
        buildable: false
        externals:
        - spec: intel-mpi@2019.0.9
          modules:
          - impi/19.0.9
      all:
        providers:
          mpi: [intel-mpi]


Additional already installed modules can be discovered and added with:

.. code-block:: bash

    spack external find

Finally, install seissol-env with

.. code-block:: bash

    spack install -j 16 seissol-env %intel ^intel-mpi

We then create the required modules with:

.. code-block:: bash

    spack module tcl refresh $(spack find -d --format "{name}{/hash:5}" seissol-env)

To access the module at start up, add to your ``~/.bashrc``:

.. code-block:: bash

    module use $SPACK_ROOT/share/spack/modules/linux-centos7-cascadelake/

Finally, install SeisSol with cmake, as usual, with ``-DHOST_ARCH=skx``.

For large runs, it is recommanded to have executable, dynamic libraries, setup and outputs on scratch.
That is, $WORK is not built for large, intensive IOs, and loading the shared libraries from 8000+ nodes at once is quite intensive.
This could potentially break the filesystem.
The dynamic libraries can be copied to $SCRATCH with the following commands:

.. code-block:: bash

    # replace by the path to your seissol executable
    mkdir -p $SCRATCH/libdump  && ldd SeisSol_Release_dskx_6_elastic | grep -E "/work|/scratch" | awk '{print $(NF-1)}' | xargs -I _ cp _ $SCRATCH/libdump

Then you can unload the seissol-env module and add the required dynamic libraries, e.g. with:

.. code-block:: bash

    export LD_LIBRARY_PATH=$SCRATCH/libdump/:$LD_LIBRARY_PATH
    module unload seissol-env

Finally, we provide an example of launch script used for running a full-machine frontera run.
In particular, note how timeout and retry count are increased.

.. code-block:: bash

    #!/bin/bash
    #SBATCH --chdir=./
    #SBATCH -o ./%j.out       # Name of stdout output file
    #SBATCH -e ./%j.out       # Name of stderr error file
    #SBATCH -p debug         # Queue (partition) name
    #SBATCH --nodes=8192
    #SBATCH --ntasks-per-node=2
    #SBATCH -t 24:00:00        # Run time (hh:mm:ss)
    #SBATCH -A EAR22007       # Project/Allocation name (req'd if you have more than 1)

    # Any other commands must follow all #SBATCH directives...
    module list
    pwd
    date

    #Prevents errors such as experience in Issue #691
    export I_MPI_SHM_HEAP_VSIZE=32768

    export OMP_NUM_THREADS=27
    export OMP_PLACES="cores(27)"
    export OMP_PROC_BIND="close"

    export XDMFWRITER_ALIGNMENT=8388608
    export XDMFWRITER_BLOCK_SIZE=8388608
    export ASYNC_MODE=THREAD
    export ASYNC_BUFFER_ALIGNMENT=8388608

    echo 'num_nodes:' $SLURM_JOB_NUM_NODES 'ntasks:' $SLURM_NTASKS
    ulimit -Ss 2097152

    source ~cazes/texascale_settings.sh
    export UCX_TLS=knem,dc
    export UCX_DC_MLX5_TIMEOUT=35000000.00us
    export UCX_DC_MLX5_RNR_TIMEOUT=35000000.00us
    export UCX_DC_MLX5_RETRY_COUNT=180
    export UCX_DC_MLX5_RNR_RETRY_COUNT=180
    export UCX_RC_MLX5_TIMEOUT=35000000.00us
    export UCX_RC_MLX5_RNR_TIMEOUT=35000000.00us
    export UCX_RC_MLX5_RETRY_COUNT=180
    export UCX_RC_MLX5_RNR_RETRY_COUNT=180
    export UCX_UD_MLX5_TIMEOUT=35000000.00us
    export UCX_UD_MLX5_RETRY_COUNT=180


    # Launch MPI code...
    seissol_exe=SeisSol_Release_dskx_6_viscoelastic2
    echo $seissol_exe
    time -p ibrun $seissol_exe parameters.par
