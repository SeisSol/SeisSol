Grace Hopper Superchip
======================

Here, we consider not necessarily a supercomputer, but the "Nvidia Grace Hopper Superchip" building block of it.

Thus, one node consists of:

- 1× Nvidia Grace CPU.
- 1× Nvidia H100 Hopper GPU.

Thus, we have the settings:

- ``HOST_ARCH=sve128``
- ``DEVICE_ARCH=sm_90``
- ``DEVICE_BACKEND=cuda``

Furthermore, we will use NVHPC where possible.

Building the Modules
~~~~~~~~~~~~~~~~~~~~

For the following, we assume that NVHPC is your default compiler. Also, we assume that you have a GCC installed as well.
To begin with, let us set the default compilers, if not done so already.

.. code-block:: bash
    export CC=nvc
    export CXX=nvc++

We will also include the build instructions for Hdf5 and Netcdf here. However, note that these are usually already pre-built on your cluster—in which case we advise you to use the available modules.

Next, we set up Python.

.. code-block:: bash

    python -m venv $SEISSOL_PREFIX
    source $SEISSOL_PREFIX/bin/activate
    pip install git+https://github.com/SeisSol/PSpaMM.git
    pip install git+https://github.com/SeisSol/gemmforge.git
    pip install git+https://github.com/SeisSol/chainforge.git

For all other modules, the default instructions suffice.

.. code-block:: bash

    wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-16.0.6.tar.gz
    tar -xf llvmorg-16.0.6.tar.gz
    mkdir -p llvm-project-llvmorg-16.0.6/build
    cd llvm-project-llvmorg-16.0.6/build
    CC=$(which gcc) CXX=$(which g++) cmake ../llvm -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;compiler-rt;openmp;polly" -DGCC_INSTALL_PREFIX="${GCC_ROOT}" -DCUDA_TOOLKIT_ROOT_DIR="${CUDA_ROOT}" -DLLVM_TARGETS_TO_BUILD="AArch64;NVPTX" -GNinja
    ninja -j $NPROCS install
    cd ../..

    wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz
    tar -xf boost_1_80_0.tar.gz
    cd boost_1_80_0

    CC=$(which gcc) CXX=$(which g++) ./bootstrap.sh --prefix=$SEISSOL_PREFIX --with-toolset=gcc --with-libraries=fiber,context,atomic,filesystem --show-libraries

    CC=$(which gcc) CXX=$(which g++) ./b2 install toolset=gcc threading=multi variant=release link=shared visibility=hidden --with-fiber --with-context --with-atomic --with-filesystem --prefix=$SEISSOL_PREFIX

.. code-block:: bash

    git clone --branch v24.02.0 --depth 1 https://github.com/AdaptiveCpp/AdaptiveCpp.git
    mkdir -p AdaptiveCpp/build
    cd AdaptiveCpp/build
    CC=$(which gcc) CXX=$(which g++) cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DWITH_OPENCL_BACKEND=OFF -DWITH_ROCM_BACKEND=OFF -DWITH_SSCP_COMPILER=OFF -DWITH_STDPAR_COMPILER=OFF -DWITH_ACCELERATED_CPU=OFF -DWITH_CUDA_BACKEND=ON -DDEFAULT_TARGETS=cuda:sm_90 -GNinja -DWITH_CUDA_NVCXX_ONLY=ON -DWITH_LEVEL_ZERO_BACKEND=OFF
    ninja -j $NPROCS install
    cd ../..

** NOTE: currently, AdaptiveCpp will not work with NVHPC out of the box for SeisSol. The reason for that is a long-standing bug in NVHPC. **


.. code-block:: bash

    wget https://ftp.mcs.anl.gov/pub/pdetools/spack-pkgs/parmetis-4.0.3.tar.gz
    tar -xvf parmetis-4.0.3.tar.gz
    cd parmetis-4.0.3
    sed -i 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g'  ./metis/include/metis.h
    make config cc=mpicc cxx=mpicxx prefix=$SEISSOL_PREFIX
    make -j $NPROCS install
    cp build/Linux-aarch64/libmetis/libmetis.a $SEISSOL_PREFIX/lib
    cp metis/include/metis.h $SEISSOL_PREFIX/include
    cd ..

    wget https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz
    tar -xf 0.8.0.tar.gz
    mkdir -p yaml-cpp-0.8.0/build
    cd yaml-cpp-0.8.0/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja
    ninja -j $NPROCS install
    cd ../..

    git clone --recursive --depth 1 https://github.com/TUM-I5/ASAGI
    mkdir -p ASAGI/build
    cd ASAGI/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja
    ninja -j $NPROCS install
    cd ../..

    wget https://www.lua.org/ftp/lua-5.4.6.tar.gz
    tar -xf lua-5.4.6.tar.gz
    cd lua-5.4.6
    make -j $NPROCS all install INSTALL_TOP=$SEISSOL_PREFIX
    cd ..

    git clone --recursive --depth 1 https://github.com/seissol/easi
    mkdir -p easi/build
    cd easi/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -DASAGI=ON -DLUA=ON -DIMPALAJIT=OFF -DEASICUBE=OFF -GNinja
    ninja -j $NPROCS install
    cd ../..

    wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    tar -xf eigen-3.4.0.tar.gz
    mkdir -p eigen-3.4.0/build
    cd eigen-3.4.0/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -GNinja
    ninja -j $NPROCS install
    cd ../..

    git clone --depth 1 https://github.com/hfp/libxsmm
    cd libxsmm
    ninja $NPROCS install
    cd ..

Compiling SeisSol
~~~~~~~~~~~~~~~~~

Finally, it's time to clone SeisSol and build it.

In total, we get the following:

.. code-block:: bash

    git clone --recursive https://github.com/SeisSol/SeisSol.git seissol
    mkdir -p seissol/build
    cd seissol/build
    cmake .. -GNinja -DPRECISION=single -DDEVICE_BACKEND=cuda -DDEVICE_ARCH=sm_90 -DHOST_ARCH=sve128 -DORDER=4 -DASAGI=ON -DNUMA_AWARE_PINNING=ON -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
    ninja

Optionally, you can install SeisSol to ``$SEISSOL_PREFIX``.
