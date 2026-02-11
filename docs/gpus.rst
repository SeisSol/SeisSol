..
  SPDX-FileCopyrightText: 2021 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

SeisSol with GPUs
=======================================


General
~~~~~~~

.. _gpu_process_pinning:

The current GPU version of SeisSol targets the latest NVIDIA, AMD, and Intel graphics cards. Therefore, you
need to have at least CUDA 11.8 or ROCm 6.0 installed in your environment. Moreover, for the best performance, make sure
that your installed MPI implementation is GPU-aware. If you are unsure, consult the documentation of the cluster you are using.

You can also utilize `UCX communication layer
<https://github.com/openucx/ucx>`_ if you need to manually configure CUDA-Aware MPI for a cluster or your local server
with an open-source MPI implementation e.g., OpenMPI.

GPU version of SeisSol follows *single rank/single GPU* strategy. Therefore,
if you want to run SeisSol on **M** nodes where each node is equipped with **N** GPUs then
make sure that you launch SeisSol with **M x N** MPI processes.

To achieve the most efficient CPU-to-GPU communication and vice versa you have
to pin your MPI processes to CPU cores which are the closest to the target
GPUs. This problem is also known as GPU affinity. Latest versions of workload
managers (e.g. SLURM) are aware of this problem and try to provide an
automatic, GPU-aware process pinning. Consider the following SLURM options:

- ``--ntasks-per-gpu``
- ``--gpu-bind``

You can also enforce good GPU affinity with rankfiles if your GPU cluster or local server
does not use a workload manager but is equipped with multiple GPUs per node.

Supported SeisSol features
~~~~~~~~~~~~~~~~~~~~~~~~~~

The SeisSol GPU version supports everything the SeisSol CPU version supports,
except for poroelasticity (which is currently being ported).

In some cases, the features are still considered "beta" due to limited testing so far;
see the following list:

- acoustic, elastic (isotropic, anisotropic), and visco-elastic wave propagation models (all stable)
- elastic-acoustic interaction (stable)
- fused simulations (beta)
- kinematic point sources (stable)
- dynamic rupture (stable; except TP in beta)
- off-fault plasticity (stable)

Compilation
~~~~~~~~~~~

To start off, make sure that you already have GPU code generator installed on your system.
If you don't have then follow this :ref:`link <gemmforge_installation>`.

After that, get the latest version of SeisSol

.. code-block:: bash

   git clone --recurse-submodules https://github.com/SeisSol/SeisSol.git seissol-gpu

Compile SeisSol with (e.g.)

.. code-block:: bash

    mkdir -p seissol-gpu/build && cd seissol-gpu/build
    cmake -DDEVICE_BACKEND=cuda -DDEVICE_ARCH=sm_70 -DHOST_ARCH=skx \
    -DCMAKE_BUILD_TYPE=Release -DPRECISION=double ..
    make -j

The following two CMake options can be useful to improve performance:

* ``USE_GRAPH_CAPTURING``: enables CUDA, HIP, or SYCL (oneAPI) graphs. These are used to speed up the kernel execution for wave propagation equations.
* ``PREMULTIPLY_FLUX``: enables the pre-multiplying of flux matrices (it was disabled for CPUs to free up cache space). This usually results in a speedup for AMD and Nvidia GPUs. By default, it is switched on when compiling for an AMD or Nvidia GPU and switched off in all other cases.
* ``DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS``: enables a hand-written kernel to speed up some internal, heavily memory-bound computations. Enabled for AMD and NVIDIA GPUs by default; but it works on all others as well.

.. _gpu-env:

Execution
~~~~~~~~~

The launching process of the GPU version of SeisSol is similar as the one of the CPU version.

.. code-block:: bash

    mpirun -n <M x N> ./SeisSol_dsm70_cuda_* ./parameters.par

The following device-specific environment variable is supported right now:

- ``SEISSOL_USM``

- ``SEISSOL_USM_MPI``

- ``SEISSOL_L2_COMPRESS``

- ``SEISSOL_TRANSFER_MODE``

``SEISSOL_USM`` specifies if the data buffers are allocated using unified/managed (i.e. CPU-accessible) memory,
or GPU memory. It is on by default on systems like the Grace Hopper Superchip or APUs like the MI300A,
and disabled on all other systems (see TODO for more information).

``SEISSOL_USM_MPI`` specifies the allocation mode for the buffers that are also used for MPI transfer.
E.g. some MPI implementations, even if GPU-aware, will treat unified/managed memory buffers are CPU buffers
otherwise.

``SEISSOL_L2_COMPRESS`` enables L2 compression on those GPUs and frameworks where it is available.
Currently, that is only CUDA; we recommend to check the respective guides for availability
(generally, Hopper or newer should always support it; and Ampere slightly restricted).
Disabled by default, and currently requires ``SEISSOL_USM=0`` and ``SEISSOL_USM_MPI=0``.

``SEISSOL_TRANSFER_MODE`` specifies how to copy GPU buffers via MPI.
The default value is ``direct`` which copies the data out of the GPU buffers directly.
In contrast, the ``host`` value means that the data will be copied to/from the host memory
before/after each send/receive operation.
That is especially useful, since some MPI implementations are not GPU-aware and do not support direct point-to-point
communication on device buffers.
As a (subpar) alternative, you can also try using ``SEISSOL_USM_MPI=1`` and ``direct`` to utilize unified/managed memory.

.. figure:: figures/gpu-comm-layer-data-flow.png
   :alt: Data Flow Diagram
   :width: 10.0cm
   :align: center

Alternative values for ``SEISSOL_TRANSFER_MODE`` are given by ``ccl`` and ``shmem``
which use NCCL/RCCL/oneCCL and NVSHMEM/ROCSHMEM/ISHMEM, respectively.
However, note that both of these are considered experimental.
