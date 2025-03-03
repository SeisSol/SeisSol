..
  SPDX-FileCopyrightText: 2023-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_heisenbug:

Heisenbug
=========

`heisenbug <https://www.geophysik.uni-muenchen.de/research/geocomputing/heisenbug>`_
is a computing cluster of the computational seismology group at LMU.
It is an AMD EPYC based machine with 128 cores that can run 256 threads (near) simultaneously.
It also has 2 GPGPUs (NVIDIA GeForce RTX 3090), that can be used to run the GPU version of SeisSol.
The RTX 3090 belongs to a consumer kind of graphics cards and thus does not perform well with double precision.
Therefore, it is preferable to compile SeisSol with single precision.

CPU and GPU seissol modules are available on heisenbug. They also integrate all libraries relevant for manually compiling another version of SeisSol.
It can be discovered at startup after adding the following to ``~/.bashrc``:

.. code-block:: bash

    module use /import/exception-dump/ulrich/spack/modules/linux-debian12-zen2

E.g. the GPU module can be loaded with:

.. code-block:: bash

    # load the (first in the list) seissol module compiled with cuda support
    module load $(module avail seissol/*-cuda-* | awk '/seissol/ {print $1}')

These modules have been compiled based on the develop branch of spack with the command:

.. code-block:: bash

    spack compiler find
    # deactivate conda module (if any)
    conda deactivate
    # python has an externally-managed-environment on heisenbug
    # We therefore need to create an environment, e.g. with:
    python3 -m venv /import/exception-dump/ulrich/python-venv
    export PATH=/import/exception-dump/ulrich/python-venv/bin:$PATH
    spack external find python
    # CPU version of seissol
    spack install -j 40 seissol@master convergence_order=4 dr_quad_rule=dunavant equations=elastic precision=single %gcc@12 ^easi +python
    # GPU version of seissol
    spack install -j 40 seissol@master convergence_order=4 dr_quad_rule=dunavant equations=elastic precision=single %gcc@12 +cuda cuda_arch=86  ^easi +python
    spack module tcl refresh $(spack find -d --format "{name}{/hash:5}" seissol)


As there is no queuing system on heisenbug, you need to make sure that nobody is running anything on the GPUs.
You can check that by running ``nvidia-smi`` (it should return ``No running processes found``).

To run on one GPU (here with order 4, elastic), use:

.. code-block:: bash

    export OMP_STACKSIZE=16M
    export MP_SINGLE_THREAD=no
    unset KMP_AFFINITY
    export OMP_PLACES="cores"
    export OMP_PROC_BIND=spread
    export XDMFWRITER_ALIGNMENT=4096
    export XDMFWRITER_BLOCK_SIZE=4096
    export ASYNC_BUFFER_ALIGNMENT=4096
    source /etc/profile.d/modules.sh
    export SEISSOL_ASAGI_MPI_MODE=OFF
    ulimit -Ss 2097152
    export OMP_NUM_THREADS=10
    export OMP_PLACES="cores"
    export OMP_PROC_BIND=spread
    seissol-launch SeisSol_Release_ssm_86_cuda_4_elastic ./parameters.par

(`seissol-launch` is a simple bash helper script).

On 2 ranks, use:

.. code-block:: bash

    export OMP_STACKSIZE=16M
    export MP_SINGLE_THREAD=no
    unset KMP_AFFINITY
    export OMP_PLACES="cores"
    export OMP_PROC_BIND=spread
    export XDMFWRITER_ALIGNMENT=4096
    export XDMFWRITER_BLOCK_SIZE=4096
    export ASYNC_BUFFER_ALIGNMENT=4096
    source /etc/profile.d/modules.sh
    export SEISSOL_ASAGI_MPI_MODE=OFF
    ulimit -Ss 2097152
    # Note that it is possible to increase OMP_NUM_THREADS
    # This will speed up (the rare) portions of the code running only CPUs, e.g. the wiggle factor calculation
    export OMP_NUM_THREADS=10
    export OMP_PLACES="cores"
    export OMP_PROC_BIND=spread
    mpirun -n 2 --map-by ppr:2:numa:pe=$OMP_NUM_THREADS --report-bindings seissol-launch SeisSol_Release_ssm_86_cuda_4_elastic ./parameters.par
