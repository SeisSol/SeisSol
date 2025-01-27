..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Frontier
========

[NOTE: this is almost a copy of the LUMI page]

Website: https://www.olcf.ornl.gov/frontier/

Guide: https://docs.olcf.ornl.gov/systems/frontier_user_guide.html

Here, we concern ourselves with running SeisSol on Frontier.

Each node consists of:

- 1× AMD Epyc 7A53 (Zen 3) CPU, configured with 4 NUMA domains
- 4× AMD Instinct MI250x GPUs, thus 8 GCDs in total

Due to the 8 GCDs, we will launch SeisSol with 8 processes per node. The architecture settings we will need for SeisSol are
``milan`` for the CPU architecture (optimizing for Zen 3), and ``gfx90a`` for the GPU architecture (targeting the MI250X).
As device backend, we use HIP, and for the SYCL implementation, we use AdaptiveCpp.

Installing Modules (without Spack)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here, we go for a build using amdclang and AdaptiveCpp. We begin by setting up an environment. Firstly, choose a folder you want to install SeisSol to and navigate to it.
Run ``pwd`` and copy the path there. Run the following script there.

.. code-block:: bash

    export SEISSOL_BASE=$(pwd)
    export SEISSOL_PREFIX=$SEISSOL_BASE/local
    export CMAKE_PREFIX_PATH=$SEISSOL_PREFIX:$CMAKE_PREFIX_PATH
    export LD_LIBRARY_PATH=$SEISSOL_PREFIX/lib:$SEISSOL_PREFIX/lib64:$LD_LIBRARY_PATH

    mkdir -p $SEISSOL_PREFIX

Next, we load the necessary modules for our SeisSol build.
We set the compilers to the cray compiler wrappers (which in our case use ``amdclang`` internally).

.. code-block:: bash

    module load Core/24.07

    module load craype-x86-trento
    module load craype-accel-amd-gfx90a

    module load PrgEnv-amd/8.5.0
    module switch amd/6.0.0
    module load rocm/6.0.0

    module load cray-hdf5-parallel/1.12.2.9
    module load cray-netcdf-hdf5parallel/4.9.0.9
    module load cray-python/3.11.5

    module load ninja
    module load cmake

    export CC=cc
    export CXX=CC
    export FC=ftn

    export HIP_PATH=/opt/rocm-6.0.0

We also require a small hotfix for pkg-config, as required by easi (and subsequently also SeisSol) right now. It is needed to work correctly with the Cray environment (only the folders ``hdf5-parallel`` and ``netcdf-hdf5parallel`` are included by default; but these do not contain the non-parallel pkgconfigs):

.. code-block:: bash

    export PKG_CONFIG_PATH=/opt/cray/pe/netcdf/4.9.0.9/amd/5.0/lib/pkgconfig:/opt/cray/pe/hdf5/1.12.2.7/amd/5.0/lib/pkgconfig:$PKG_CONFIG_PATH

Next, we also start up our Python installation. The virtual environment sets additional paths for e.g. executables to our prefix directory automatically.

.. code-block:: bash

    python -m venv $SEISSOL_PREFIX
    source $SEISSOL_PREFIX/bin/activate
    pip install setuptools
    pip install numpy
    pip install git+https://github.com/SeisSol/PSpaMM.git
    pip install git+https://github.com/SeisSol/gemmforge.git
    pip install git+https://github.com/SeisSol/chainforge.git

Then, we can start installing the modules. For AdaptiveCpp, we need a suitable Boost installation. That can be accomplished as follows:

.. code-block:: bash

    wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz
    tar -xf boost_1_80_0.tar.gz
    cd boost_1_80_0

    ./bootstrap.sh --prefix=$SEISSOL_PREFIX --with-toolset=gcc --with-libraries=fiber,context,atomic,filesystem --show-libraries

    ./b2 install toolset=gcc threading=multi variant=release link=shared visibility=hidden --with-fiber --with-context --with-atomic --with-filesystem --prefix=$SEISSOL_PREFIX

    cd ..

Next, we build AdaptiveCpp. Note that we need to switch off everything but ROCm for the installation to work smoothly.

.. code-block:: bash

    git clone --branch v23.10.0 --depth 1 https://github.com/AdaptiveCpp/AdaptiveCpp.git
    mkdir -p AdaptiveCpp/build
    cd AdaptiveCpp/build
    cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DWITH_OPENCL_BACKEND=OFF -DWITH_ROCM_BACKEND=ON -DWITH_SSCP_COMPILER=OFF -DWITH_STDPAR_COMPILER=OFF -DWITH_ACCELERATED_CPU=OFF -DWITH_CUDA_BACKEND=OFF -DWITH_LEVEL_ZERO_BACKEND=OFF -DDEFAULT_TARGETS=hip:gfx90a
    ninja install
    cd ../..

The rest of the packages can be installed as usual.

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

For ImpalaJIT (depending on the former two):

.. code-block:: bash

    git clone https://github.com/uphoffc/ImpalaJIT.git
    cd ImpalaJIT/
    mkdir -p build
    cd build/
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja
    ninja install
    cd ../..

For easi (depending on the former two):

.. code-block:: bash

    git clone --recursive --depth 1 https://github.com/seissol/easi
    mkdir -p easi/build
    cd easi/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja -DASAGI=ON -DLUA=ON -DIMPALAJIT=ON -DEASICUBE=OFF
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

In case there are problems with using libxsmm, you can also consider using only PSpaMM instead; at a tiny performance penalty.

Compiling SeisSol
~~~~~~~~~~~~~~~~~

Finally, it's time to clone SeisSol and build it.

However, we need to apply a small hotfix here, since the Cray compiler environment does not work with AdaptiveCpp (it causes problems with finding MPI, the filesystem headers etc.). As a workaround, we compile SeisSol with ``amdclang`` directly, and add the necessary flags from the Cray environment as compiler flags (that can be done by ``CC --cray-print-opts=all``, the same with ``cc`` and ``ftn``). There is currently a problem with Graph Capturing, so we turn it off for the time-being

In total, we get the following:

.. code-block:: bash

    git clone --recursive https://github.com/SeisSol/SeisSol.git seissol
    mkdir -p seissol/build
    cd seissol/build
    CC=amdclang CXX=amdclang++ CFLAGS=$(cc --cray-print-opts=all) CXXFLAGS=$(CC --cray-print-opts=all) cmake .. -GNinja -DPRECISION=single -DDEVICE_BACKEND=hip -DDEVICE_ARCH=gfx90a -DHOST_ARCH=milan -DORDER=4 -DASAGI=ON -DNUMA_AWARE_PINNING=ON -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DUSE_GRAPH_CAPTURING=OFF
    ninja

Optionally, you can install SeisSol to ``$SEISSOL_PREFIX``.

Running Jobs
~~~~~~~~~~~~

Attached is a job script which does the pinning for us.
The pinning on the Frontier nodes needs some special attention, since 8 out of the 64 cores are reserved for the OS.

.. code-block:: bash

    #!/usr/bin/env bash
    #SBATCH --job-name=seissol   # Job name
    #SBATCH --nodes=<NUMBER-OF-NODES>               # Total number of nodes
    #SBATCH --account=<your-project>  # Project for billing
    #SBATCH --mail-user=<your-mail>
    #SBATCH --time=01:00:00       # Run time (d-hh:mm:ss)
    #SBATCH --output=seissol-output.log # Name of stdout output file
    #SBATCH --error=seissol-error.log  # Name of stderr error file
    #SBATCH --partition=standard-g  # Partition (queue) name
    #SBATCH --ntasks-per-node=8     # 8 MPI ranks per node
    #SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
    #SBATCH --mail-type=all         # Send email at begin and end of job
    #SBATCH --exclusive
    #SBATCH --requeue

    cat << EOF > select_gpu
    #!/bin/bash

    export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
    exec \$*
    EOF

    chmod +x ./select_gpu

    CPU_BIND="7e000000000000,7e00000000000000"
    CPU_BIND="${CPU_BIND},7e0000,7e000000"
    CPU_BIND="${CPU_BIND},7e,7e00"
    CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

    export MPICH_GPU_SUPPORT_ENABLED=1
    export HSA_XNACK=0

    export OMP_NUM_THREADS=3
    export OMP_PLACES="cores(3)"
    export OMP_PROC_BIND=close

    export DEVICE_STACK_MEM_SIZE=4
    export SEISSOL_FREE_CPUS_MASK="52-54,60-62,20-22,28-30,4-6,12-14,36-38,44-46"

    srun --cpu-bind=mask_cpu:${CPU_BIND} ./select_gpu ./SeisSol_Release_sgfx90a_hip_6_elastic parameters.par
