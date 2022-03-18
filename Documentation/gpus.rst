SeisSol with GPUs
=======================================


General
~~~~~~~

.. _gpu_process_pinning:

The current GPU version of SeisSol targets latest NVidia Graphics cards. Therefore, you
need to have at least **CUDA 10** installed in your environment. Moreover, make sure
that your installed MPI implementation is **CUDA-Aware**. Here is a list of 
known CUDA-Aware MPI vendors:

- OpenMPI
- MVAPICH2
- CRAY MPI
- IBM Platform MPI
- SGI MPI

Please, referer to a corresponding documentation if you need to installed
CUDA-Aware MPI manually. We also ecourage you to bind your MPI with `UCX communication layer
<https://github.com/openucx/ucx>`_ if you need to manually configure CUDA-Aware MPI for a cluster or your local server 
with an open-source MPI implementation e.g., OpenMPI.

GPU version of SeisSol follows *single rank/single GPU* strategy. Therefore, 
if you want to run SeisSol on **M** nodes where each node is equipped with **N** GPUs then
make sure that you launch SeisSol with **M x N** MPI processes. 

To achieve the most efficient CPU-to-GPU communication and vice versa you have 
to pin your MPI processes to CPU cores which are the closest to the target 
GPUs. This problem is also known as GPU affinity. Latest versions of workload 
managers (e.g., SLURM) are aware about this problem and try to provide an 
automatic, GPU-aware process pinning. Consider the following SLURM options:

- `--ntasks-per-gpu`
- `--gpu-bind`

You can also enforce good GPU affinity with rankfiles if your GPU cluster or local server
does not use a workload manager but is equipped with multiple GPUs per node.

.. figure:: LatexFigures/GpuCpuProcessPinning.png
   :alt: Process Pinning
   :width: 16.00000cm
   :align: center

   Correct process pinning of 4 MPI processes where each process
   controls 3 OpenMP threads and one communication thread.

Supported SeisSol Features
~~~~~~~~~~~~~~~~~~~~~~~~~~

- elastic wave propagation model with kinematic point sources
- off-fault plasticity model


Compilation
~~~~~~~~~~~

To start off, make sure that you already have **GemmForge** installed on your system. 
If you don't have then follow this :ref:`link <gemmforge_installation>`.

After that, get the latest version of SeisSol

.. code-block:: bash

   git clone --recurse-submodules https://github.com/SeisSol/SeisSol.git seissol-gpu

Compile SeisSol with (e.g.)

.. code-block:: bash

    mkdir -p seissol-gpu/build && cd seissol-gpu/build 
    cmake -DDEVICE_BACKEND=cuda -DDEVICE_ARCH=sm_70 -DHOST_ARCH=skx \
    -DCOMMTHREAD=ON -DCMAKE_BUILD_TYPE=Release -DPRECISION=double ..
    make -j

Execution
~~~~~~~~~

A launching process of GPU version of SeisSol is exactly the same as in case
of the CPU version. 

.. code-block:: bash

    mpirun -n <M x N> ./SeisSol_dsm70_cuda_* ./parameters.par

It is important to know that the GPU version of SeisSol allocates 1GB of
GPU memory at the beginning of SeisSol execution. It is necessary for fast allocation/deallocation
of GPU memory needed for holding temporary data. The default value can be changed by setting
a necessary one to **DEVICE_STACK_MEM_SIZE** environment variable. For example,
the following will force SeisSol to allocate 1.5GB of stack GPU memory for temporary data:


.. code-block:: bash
    
    export DEVICE_STACK_MEM_SIZE=1.5
    mpirun -n <M x N> ./SeisSol_dsm70_cuda_* ./parameters.par


heisenbug
~~~~~~~~~

`heisenbug <https://www.geophysik.uni-muenchen.de/research/geocomputing/heisenbug>`_ is a computing cluster of the computational seismology group at LMU.
It is an AMD EPYC based machine with 128 cores that can run 256 threads (near) simultaneously. 
It also has 2 GPGPUs (NVIDIA GeForce RTX 3090), that can be used to run the GPU version of SeisSol.


A CUDA-Aware MPI on heisenbug is loaded with:

.. code-block:: bash

    module load gcc/10.2.0
    module load mpi.ompi/4.1.0

GemmForge and ChainForge (for GPU) can be installed with:

.. code-block:: bash

    pip3 install git+https://github.com/ravil-mobile/gemmforge.git
    pip3 install git+https://github.com/ravil-mobile/chainforge.git

All dependencies required for the CPU version of SeisSol should also be installed.
A folder containing all precompiled requirements is available at:

.. code-block:: bash

    /export/dump/ulrich/SeisSol-GPU-build

To compile the GPU version of SeisSol on heisenbug, use the following cmake options ``-DDEVICE_ARCH=sm_86 -DHOST_ARCH=hsw``.
Use ``-DCOMMTHREAD=ON`` for multiple GPUs, and ``-DCOMMTHREAD=OFF`` for one GPU.

To run on one GPU (here with order 4, elastic), use simply:

.. code-block:: bash

    export OMPI_COMM_WORLD_LOCAL_RANK=0
    export OMPI_COMM_WORLD_LOCAL_SIZE=1
    export OMPI_MCA_btl=vader,self
    SeisSol_Release_ssm_86_cuda_4_elastic parameters.par

On 2 ranks, use:

.. code-block:: bash

    export OMPI_MCA_btl=vader,self
    mpirun -n 2 SeisSol_Release_ssm_86_cuda_4_elastic parameters.para

