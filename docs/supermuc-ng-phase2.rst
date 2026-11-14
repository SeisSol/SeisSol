..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

SuperMUC NG Phase 2
===================

Website: https://doku.lrz.de/supermuc-ng-10745965.html

The nodes consist of (cf. https://doku.lrz.de/hardware-of-supermuc-ng-phase-2-222891050.html ):

- 2× Intel Xeon Platinum 8480L CPUs (Sapphire Rapids, 56 cores each), 512 GB of memory
- 4× Intel Data Center GPU Max 1550 GPUs (Ponte Vecchio), with two stacks each.

With ``ZE_FLAT_DEVICE_HIERARCHY=FLAT``, each stack is a GPU of its own: a node then has eight,
and we run one MPI rank per stack.
We use ``skx`` as CPU architecture, ``pvc`` as GPU architecture, and ``oneapi`` as backend (and compiler).

These instructions were tested in November 2026. The GPU driver of the nodes at that time
(compute-runtime 23.43 with the Level Zero loader 1.14) limits the oneAPI versions that work, see
:ref:`supermuc-ng-phase2-known-issues`: use oneAPI 2025.0.

The login nodes have no access to the internet. Download the sources (and Python packages) elsewhere and copy them over,
or see :ref:`pypi_behind_firewall`.

Environment
~~~~~~~~~~~

The spack stack 26.0.0 of the LRZ has oneAPI 2025.0 and most libraries SeisSol needs.
Its HDF5 and ParMETIS are built for Intel MPI 2021.16, whose modules conflict with the Intel MPI of oneAPI 2025.0
(2021.14, the same ABI), so we take them by their paths.
Put the following into a file, e.g. ``env.sh``, and source it for building and in the job scripts:

.. code-block:: bash

    module switch stack/26.0.0
    module load intel-toolkit/2025.0.1
    module load numactl/2.0.19-intel25 metis/5.1.0-intel25-i64-r64 eigen/5.0.1-intel25
    module load cmake/3.31.11 ninja/1.13.2 libxsmm/1.17-intel25-shared python/3.14.0

    # the install prefix of a module of the stack
    prefix() { module show "$1" 2>&1 | grep -oE "_BASE [^ ]+" | head -1 | cut -d" " -f2; }
    export HDF5_ROOT=$(prefix hdf5/1.14.6-intel25-impi)
    export PARMETIS_DIR=$(prefix parmetis/4.0.3-intel25-impi-i64-r64)
    export METIS_DIR=$METIS_BASE

    export SEISSOL_PREFIX=<your-install-folder>
    export PATH=$SEISSOL_PREFIX/bin:$HDF5_ROOT/bin:$PATH
    export LD_LIBRARY_PATH=$SEISSOL_PREFIX/lib:$SEISSOL_PREFIX/lib64:$HDF5_ROOT/lib:$PARMETIS_DIR/lib:$LD_LIBRARY_PATH
    export CMAKE_PREFIX_PATH=$SEISSOL_PREFIX:$HDF5_ROOT:$PARMETIS_DIR:$METIS_DIR:$CMAKE_PREFIX_PATH
    export CC=icx CXX=icpx

The Python 3.14 of the LRZ comes with NumPy and PyYAML, which the code generators need. We make a virtual
environment on top of it, with PSpaMM and TensorForge (the code generator for the GPU kernels) installed from their sources,
i.e. clones of https://github.com/SeisSol/PSpaMM and https://github.com/SeisSol/TensorForge copied to the machine.
(The Tiny Tensor Compiler, the other code generator for Intel GPUs, see :ref:`build_dependencies`, we have not tried there.)

.. code-block:: bash

    python3 -m venv --system-site-packages $SEISSOL_PREFIX
    source $SEISSOL_PREFIX/bin/activate
    pip install --no-index --no-deps --no-build-isolation ./PSpaMM ./tensorforge

(Add ``source $SEISSOL_PREFIX/bin/activate`` to ``env.sh`` as well.)

Libraries
~~~~~~~~~

Eigen, METIS/ParMETIS, HDF5 and the libxsmm generator come from the modules above. What remains are
yaml-cpp, Lua and easi (with Lua, which also evaluates ``!FunctionMap``), from their release archives:

.. code-block:: bash

    cmake -S yaml-cpp-0.8.0 -B yaml-cpp-0.8.0/build -GNinja -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    ninja -C yaml-cpp-0.8.0/build install

    make -C lua-5.4.6 linux CC=icx MYCFLAGS=-fPIC
    make -C lua-5.4.6 install INSTALL_TOP=$SEISSOL_PREFIX

    cmake -S easi -B easi/build -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX \
      -DASAGI=OFF -DLUA=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    ninja -C easi/build install

NetCDF (for ASAGI and for point sources in the NRF format) is available as ``netcdf-c/4.9.3-intel25-impi``,
by its path like HDF5; we have not tested it there.

Compiling SeisSol
~~~~~~~~~~~~~~~~~

Compiling the GPU kernels ahead of time for the Max 1550 takes a while: build on a compute node, e.g. in a job
of one node with ``ninja -j 96``, which takes about seven minutes.

.. code-block:: bash

    git clone --recursive https://github.com/SeisSol/SeisSol.git seissol
    cmake -S seissol -B seissol/build -GNinja -DCMAKE_BUILD_TYPE=Release \
      -DDEVICE_BACKEND=oneapi -DDEVICE_ARCH=pvc -DHOST_ARCH=skx \
      -DPRECISION=single -DORDER=6 -DEQUATIONS=elastic -DNETCDF=OFF -DASAGI=OFF \
      -DUSE_GRAPH_CAPTURING=OFF \
      -DTensorForge_DIR=$(python3 -c "import tensorforge, os; print(os.path.dirname(tensorforge.__file__))")/share/cmake
    ninja -C seissol/build -j 96

The build directory then holds ``seissol-elastic-o6-f32`` and ``seissol-launch``.
``USE_GRAPH_CAPTURING=OFF`` (an option of the Device library) leaves out the SYCL graphs, see the known issues below.

Running Jobs
~~~~~~~~~~~~

A job script for eight MPI ranks per node, one per GPU stack. ``seissol-launch`` leaves each rank the stack of its
local rank (``ZE_AFFINITY_MASK``); Intel MPI pins each rank to an eighth of the node, i.e. 14 cores next to its GPU.

.. code-block:: bash

    #!/bin/bash
    #SBATCH -J seissol
    #SBATCH -A <your-project>
    #SBATCH -p general
    #SBATCH -N <number-of-nodes>
    #SBATCH --ntasks-per-node=8
    #SBATCH -t 01:00:00
    #SBATCH --export=NONE
    #SBATCH --get-user-env

    module load slurm_setup
    source env.sh

    # (LRZ's mpi_settings select the psm3 provider with GPU direct, which the halo exchange between nodes needs)
    export ZE_FLAT_DEVICE_HIERARCHY=FLAT

    export I_MPI_PIN_DOMAIN=auto
    export OMP_NUM_THREADS=12 OMP_PLACES=cores OMP_PROC_BIND=close

    # the host actor, without a communication thread (see the time stepping on GPUs)
    export SEISSOL_HOST_ACTOR=1 SEISSOL_COMMTHREAD=0

    mpiexec -n $SLURM_NTASKS -ppn 8 ./seissol-launch ./seissol-elastic-o6-f32 parameters.par

The host actor (:ref:`gpu-time-stepping`) made the runs on the Max 1550 faster; concurrent clusters
(``SEISSOL_CONCURRENT_CLUSTERS=1``) did not change much there, and the time-stepping plan made them slower.
Built as above, the verification cases of single precision pass, and two runs of the same binary came out bitwise
the same, also with these options.

.. _supermuc-ng-phase2-known-issues:

Known issues
~~~~~~~~~~~~

As of November 2026, with the GPU driver of the nodes (compute-runtime 23.43, Level Zero loader 1.14):

- oneAPI 2025.3 and later find no GPU through their default Level Zero adapter (``level_zero backend failed with
  error: 37 (UR_RESULT_ERROR_UNINITIALIZED)``), only with ``SYCL_UR_USE_LEVEL_ZERO_V2=1``. On this driver, that adapter
  (in oneAPI 2025.2 and 2025.3) reports a copy as complete before the work submitted before it on the same in-order queue
  is, and 2025.2 also ran host tasks too early: the samples of the receivers then differ from run to run, and the halo
  exchange delivers stale data. With oneAPI 2026.1, a small test program crashed with it.
- With the default adapter of oneAPI 2025.2, SeisSol hangs after a few steps; that version also lists each GPU stack
  twice (16 GPUs per node in ``sycl-ls``).
- With oneAPI 2025.2, the ``psm3`` provider with GPU direct, which LRZ's ``mpi_settings`` select, crashes in
  ``MPI_Init``. ``PSM3_GPUDIRECT=0`` avoids the crash, but then the exchange of GPU memory between nodes fails
  (``reg_mr failed ... Bad address``).

- With the SYCL graphs (on by default), the results of the clusters one after another differed slightly from run to
  run, and the runs were slower; hence ``USE_GRAPH_CAPTURING=OFF`` above.

oneAPI 2025.0 (and 2024.2, as far as we tested) do not have the problems of the adapters and of ``psm3`` with this driver.
