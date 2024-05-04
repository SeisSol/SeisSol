LUMI
====

Website: https://www.lumi-supercomputer.eu/

Here, we concern ourselves with running SeisSol on the **LUMI-G** partition; that is, on the GPU partition.

The nodes then consist of:

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

    mkdir -p $SEISSOL_PREFIX

Next, we load the necessary modules for our SeisSol build.
We set the compilers to the cray compiler wrappers (which in our case use ``amdclang`` internally).

.. code-block:: bash

    module load LUMI/23.09 partition/G
    module load cpeAMD/23.09
    module load rocm/5.6.1
    module load amd/5.6.1

    module load Boost/1.82.0-cpeAMD-23.09
    module load Eigen/3.4.0

    module load cray-hdf5-parallel/1.12.2.7
    module load cray-netcdf-hdf5parallel/4.9.0.7
    module load cray-python/3.10.10

    export CC=cc
    export CXX=CC
    export FC=ftn

Next, we also start up our Python installation.

.. code-block:: bash

    python -m venv $SEISSOL_PREFIX
    source $SEISSOL_PREFIX/bin/activate
    pip install setuptools
    pip install numpy
    pip install git+https://github.com/SeisSol/PSpaMM.git
    pip install git+https://github.com/SeisSol/gemmforge.git
    pip install git+https://github.com/SeisSol/chainforge.git

Then, we can start installing the modules. For convenience, we also add Ninja as a build tool here first.

.. code-block:: bash

    git clone --branch v1.12.0 --depth 1 https://github.com/ninja-build/ninja.git
    mkdir -p ninja/build
    cd ninja/build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
    make -j 10 install
    cd ../..

Next, we choose AdaptiveCpp. Note that we need to switch off everything but ROCm for the installation to work smoothly.

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

For easi (depending on the former two):

.. code-block:: bash

    git clone --recursive --depth 1 https://github.com/seissol/easi
    mkdir -p easi/build
    cd easi/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja -DASAGI=ON -DLUA=ON -DIMPALAJIT=OFF -DEASICUBE=OFF
    ninja install
    cd ../..

For libxsmm (note that we need 1.17 sharp; the latest main will not work as intended with the generator):

.. code-block:: bash

    git clone --branch 1.17 --depth 1 https://github.com/hfp/libxsmm
    cd libxsmm
    make generator
    cp bin/libxsmm_gemm_generator $SEISSOL_PREFIX/bin
    cd ..

Compiling SeisSol
~~~~~~~~~~~~~~~~~

Finally, it's time to clone SeisSol and build it.

.. code-block:: bash

    git clone --recursive https://github.com/SeisSol/SeisSol.git seissol
    mkdir -p seissol/build
    cd seissol/build
    cmake .. -GNinja -DPRECISION=single -DDEVICE_BACKEND=hip -DDEVICE_ARCH=gfx90a -DHOST_ARCH=milan -DORDER=4 -DASAGI=ON -DNUMA_AWARE_PINNING=ON -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
    ninja

Optionally, you can install SeisSol to ``$SEISSOL_PREFIX``.

Running Jobs
~~~~~~~~~~~~

Attached is a job script which does the pinning for us.
The pinning on the LUMI nodes needs some special attention, since 8 out of the 64 cores are reserved for the OS (cf. https://lumi-supercomputer.github.io/LUMI-training-materials/User-Updates/Update-202308/lumig-lownoise/ ).

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
    export HSA_XNACK=1

    export OMP_NUM_THREADS=3
    export OMP_PLACES="cores(3)"
    export OMP_PROC_BIND=close

    export DEVICE_STACK_MEM_SIZE=4
    export SEISSOL_FREE_CPUS_MASK="52-54,60-62,20-22,28-30,4-6,12-14,36-38,44-46"

    srun --cpu-bind=mask_cpu:${CPU_BIND} ./select_gpu SeisSol_Release_sgfx90a_hip_6_elastic parameters.par
