..
  SPDX-FileCopyrightText: 2024-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Leonardo
========

Website: https://leonardo-supercomputer.cineca.eu/

We will focus here on the Booster module of Leonardo, i.e. we use GPUs. That is, one node consists of:

- 1× Intel Xeon Platinum 8358 (Sapphire Rapids) CPU @ 2 sockets × 32 cores × 2 HT
- 4× Nvidia A100 SXM6 GPUs

Thus, we will run SeisSol with 4 ranks per node. As architectures, we compile for the host/CPU architecture ``skx``, and use ``sm_80`` as architecture for the GPUs, together
with CUDA as device backend. For the SYCL parts, we use AdaptiveCpp (formerly known as hipSYCL or Open SYCL).

Build Instructions (without Spack)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here, we go for a build using gcc and AdaptiveCpp. We begin by setting up an environment. Firstly, choose a folder you want to install SeisSol to and navigate to it.
Run ``pwd`` and copy the path there. Run the following script there.

.. code-block:: bash

    export SEISSOL_BASE=$(pwd)
    export SEISSOL_PREFIX=$SEISSOL_BASE/local
    export CMAKE_PREFIX_PATH=$SEISSOL_PREFIX:$CMAKE_PREFIX_PATH

    mkdir -p $SEISSOL_PREFIX

Next, we load the necessary modules for our SeisSol build.
We load GCC explicitly, so that the environment variables for ``CC``, ``CXX``, etc. are set for us automatically.

.. code-block:: bash

    module load python/3.10.8--gcc--11.3.0
    module load cuda/11.8
    module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8
    module load gcc/11.3.0
    module load ninja
    module load hdf5/1.12.2--openmpi--4.1.4--gcc--11.3.0
    module load netcdf-c/4.9.0--openmpi--4.1.4--gcc--11.3.0

It can be useful to place these module loads into a script of their own for running jobs later-on.

Unforunately, we need to compile LLVM ourselves, as the given LLVM module does not have an NVPTX backend enabled (note: consider switching your build directory to TMP or SCRATCH, since LLVM can be big).

.. code-block:: bash

    export GCC_ROOT=$(dirname $(dirname $(which gcc)))
    export CUDA_ROOT=$(dirname $(dirname $(which nvcc)))

    wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-16.0.6.tar.gz
    tar -xf llvmorg-16.0.6.tar.gz
    mkdir -p llvm-project-llvmorg-16.0.6/build
    cd llvm-project-llvmorg-16.0.6/build
    cmake ../llvm -GNinja -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;compiler-rt;openmp;polly" -DGCC_INSTALL_PREFIX="${GCC_ROOT}" -DCUDA_TOOLKIT_ROOT_DIR="${CUDA_ROOT}" -DLLVM_TARGETS_TO_BUILD="X86;NVPTX"
    ninja install
    cd ../..

We need to install Boost manually as well, since the necessary components (fiber, context) are not included in the supplied Leonardo modules.
(note that atomic and filesystem are required for the fiber module to work—not specifying them during the b2 command will lead to the cmake files not being generated)

.. code-block:: bash

    wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz
    tar -xf boost_1_80_0.tar.gz
    cd boost_1_80_0

    ./bootstrap.sh --prefix=$SEISSOL_PREFIX --with-toolset=gcc --with-libraries=fiber,context,atomic,filesystem --show-libraries

    ./b2 install toolset=gcc threading=multi variant=release link=shared visibility=hidden --with-fiber --with-context --with-atomic --with-filesystem --prefix=$SEISSOL_PREFIX

    cd ..

For AdaptiveCpp itself, you can then do:

.. code-block:: bash

    git clone --branch v23.10.0 --depth 1 https://github.com/AdaptiveCpp/AdaptiveCpp.git
    mkdir -p AdaptiveCpp/build
    cd AdaptiveCpp/build
    cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DWITH_OPENCL_BACKEND=OFF -DWITH_ROCM_BACKEND=OFF -DWITH_SSCP_COMPILER=OFF -DWITH_STDPAR_COMPILER=OFF -DWITH_ACCELERATED_CPU=OFF -DWITH_CUDA_BACKEND=ON -DWITH_LEVEL_ZERO_BACKEND=OFF -DDEFAULT_TARGETS=cuda:sm_80
    ninja install
    cd ../..

With all these packages and AdaptiveCpp at hand, you are left to install easi (optionally with Lua and ASAGI), Eigen, as well as the code generators libxsmm, PSpaMM, gemmforge and chainforge.

METIS/ParMETIS:

.. code-block:: bash

    wget https://ftp.mcs.anl.gov/pub/pdetools/spack-pkgs/parmetis-4.0.3.tar.gz
    tar -xvf parmetis-4.0.3.tar.gz
    cd parmetis-4.0.3
    sed -i 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g'  ./metis/include/metis.h
    make config cc=mpicc cxx=mpicxx prefix=$SEISSOL_PREFIX
    make install
    cp build/Linux-x86_64/libmetis/libmetis.a $SEISSOL_PREFIX/lib
    cp metis/include/metis.h $SEISSOL_PREFIX/include
    cd ..

YAML-CPP can be installed as follows:

.. code-block:: bash

    wget https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz
    tar -xf 0.8.0.tar.gz
    mkdir -p yaml-cpp-0.8.0/build
    cd yaml-cpp-0.8.0/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja
    ninja install
    cd ../..

For easi, Eigen and libxsmm, the default instructions suffice.

For ASAGI:

.. code-block:: bash

    git clone --recursive --depth 1 https://github.com/TUM-I5/ASAGI
    mkdir -p ASAGI/build
    cd ASAGI/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja
    ninja install
    cd ../..

For LUA:

.. code-block:: bash

    wget https://www.lua.org/ftp/lua-5.4.6.tar.gz
    tar -xf lua-5.4.6.tar.gz
    cd lua-5.4.6
    make all install INSTALL_TOP=$SEISSOL_PREFIX
    cd ..

For easi (depending on the former two):

.. code-block:: bash

    git clone --recursive --depth 1 https://github.com/seissol/easi
    mkdir -p easi/build
    cd easi/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja -DASAGI=ON -DLUA=ON -DIMPALAJIT=OFF -DEASICUBE=OFF
    ninja install
    cd ../..

For Eigen:

.. code-block:: bash

    wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    tar -xf eigen-3.4.0.tar.gz
    mkdir -p eigen-3.4.0/build
    cd eigen-3.4.0/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -GNinja
    ninja install
    cd ../..

For libxsmm (note that we need 1.17 sharp; the latest main will not work as intended with the generator):

.. code-block:: bash

    git clone --branch 1.17 --depth 1 https://github.com/hfp/libxsmm
    cd libxsmm
    make generator
    cp bin/libxsmm_gemm_generator $SEISSOL_PREFIX/bin
    cd ..

For the Python modules (PSpaMM, gemmforge, chainforge), it is best to install them to a virtual environment:

.. code-block:: bash

    python -m venv $SEISSOL_PREFIX
    source $SEISSOL_PREFIX/bin/activate
    pip install git+https://github.com/SeisSol/PSpaMM.git
    pip install git+https://github.com/SeisSol/gemmforge.git
    pip install git+https://github.com/SeisSol/chainforge.git

Compiling SeisSol
~~~~~~~~~~~~~~~~~

Finally, it's time to clone SeisSol and build it.

.. code-block:: bash

    git clone --recursive https://github.com/SeisSol/SeisSol.git seissol
    mkdir -p seissol/build
    cd seissol/build
    cmake .. -GNinja -DPRECISION=single -DDEVICE_BACKEND=cuda -DDEVICE_ARCH=sm_80 -DHOST_ARCH=skx -DORDER=4 -DASAGI=ON -DNUMA_AWARE_PINNING=ON -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
    ninja

Optionally, you can install SeisSol to ``$SEISSOL_PREFIX``.

Running Jobs
~~~~~~~~~~~~

Attached is a (preliminary) job script which does the necessary pinning for us.

.. code-block:: bash

    #!/usr/bin/env bash
    #SBATCH --account=<PROJECT_NAME>
    #SBATCH --job-name=<JOB_NAME>
    #SBATCH --nodes=<NODE_COUNT>
    #SBATCH --time=<TIME>
    #SBATCH --qos=normal
    #SBATCH --partition=boost_usr_prod
    #SBATCH --ntasks-per-node=4
    #SBATCH --cpus-per-task=8
    #SBATCH --gres=gpu:4
    #SBATCH --mem=200G
    #SBATCH --exclusive
    #SBATCH --output=seissol-stdout.log
    #SBATCH --error=seissol-stderr.log
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

    SEISSOL_EXE=$(basename $(realpath ./SeisSol_Release*))

    ulimit -c unlimited

    CPU_BIND="mask_cpu"
    CPU_BIND="${CPU_BIND}:000f000f"
    CPU_BIND="${CPU_BIND},00f000f0"
    CPU_BIND="${CPU_BIND},0f000f00"
    CPU_BIND="${CPU_BIND},f000f000"

    srun --cpu-bind=${CPU_BIND} ./select_gpu ./${SEISSOL_EXE} ./parameters.par

