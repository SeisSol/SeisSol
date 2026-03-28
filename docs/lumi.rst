..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

LUMI
====

[NOTE: this is almost a copy of the Frontier page]

Website: https://www.lumi-supercomputer.eu/

Here, we concern ourselves with running SeisSol on the **LUMI-G** partition; that is, on the GPU partition.

The nodes then consist of:

- 1× AMD Epyc 7A53 (Zen 3) CPU, configured with 4 NUMA domains
- 4× AMD Instinct MI250X GPUs, thus 8 GCDs in total

Due to the 8 GCDs, we will launch SeisSol with 8 processes per node. The architecture settings we will need for SeisSol are
``milan`` for the CPU architecture (optimizing for Zen 3), and ``gfx90a`` for the GPU architecture (targeting the MI250X).
As device backend, we use HIP.

Installing Modules (without Spack)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here, we go for a build using amdclang and AdaptiveCpp. We begin by setting up an environment.
First, we create our installation root.
Please execute the code below (updating project name and user name).
We recommend using `/scratch` as the `/project` directory is typically limited in size and file count, and because the `/scratch` (to our experience) is never cleaned.

.. code-block:: bash

    # Define your base directory
    export SEISSOL_BASE=/scratch/<your_project>/<your_name>/seissol_base
    export SEISSOL_PREFIX=$SEISSOL_BASE/local

    # Create directories
    mkdir -p $SEISSOL_PREFIX/lib/pkgconfig
    # Create python environment
    python -m venv $SEISSOL_PREFIX


To allow CMake to find Cray-specific NetCDF and MPI dependencies, we create a wrapper file, netcdf.pc within our local prefix. Note that we now place it in $SEISSOL_PREFIX/lib/pkgconfig.
This hotfix is required, because the Cray packages provide only parallel pkg-config files (e.g. hdf5_parallel.pc), which list dependencies such as netcdf.
These dependency names are standard upstream pkg-config names, but Cray’s NetCDF packages do not ship matching .pc files (netcdf.pc).
Please execute the following code:

.. code-block:: bash

    cat > $SEISSOL_PREFIX/lib/pkgconfig/netcdf.pc <<'EOF'
    prefix=/opt/cray/pe/netcdf-hdf5parallel/4.9.0.17/amd/6.0
    exec_prefix=\${prefix}
    libdir=\${prefix}/lib
    includedir=\${prefix}/include

    Name: netcdf
    Description: NetCDF (parallel) wrapper
    Version: 4.9.0.17
    Requires: hdf5_hl_parallel hdf5_parallel mpich >= 7.0

    Cflags: -I\${includedir}
    Libs: -L\${libdir} -lnetcdf
    EOF


Next, we  update `~/.bashrc`, to automatically load the necessary modules for our SeisSol build, and set up installation paths.
We set the compilers to the cray compiler wrappers (which in our case use ``amdclang`` internally).
Add to your `~/.bashrc` file (again update project name and user name):

.. code-block:: bash

    module load LUMI partition/G
    module load rocm
    module load cpeAMD
    module load Eigen
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load cray-python

    export SEISSOL_BASE=/scratch/<your_project>/<your_name>/seissol_base
    export SEISSOL_PREFIX=$SEISSOL_BASE/local

    export CMAKE_PREFIX_PATH=$SEISSOL_PREFIX:$CMAKE_PREFIX_PATH
    export PKG_CONFIG_PATH=$SEISSOL_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

    if [ -f "$SEISSOL_PREFIX/bin/activate" ]; then
        source $SEISSOL_PREFIX/bin/activate
    fi

    export CC=cc
    export CXX=CC
    export FC=ftn

Don't forget to source the  `~/.bashrc` after these changes:

.. code-block:: bash

    source ~/.bashrc

Next, we also start up our Python installation. The virtual environment sets additional paths for e.g. executables to our prefix directory automatically.

.. code-block:: bash

    pip install setuptools
    pip install numpy
    pip install ninja
    pip install pybind11
    pip install git+https://github.com/SeisSol/PSpaMM.git
    pip install git+https://github.com/SeisSol/gemmforge.git
    pip install git+https://github.com/SeisSol/chainforge.git
    pip install git+https://github.com/SeisSol/TensorForge.git

Then, we can start installing the modules.
The required packages can be installed as described in the dependency installing guide itself.
For convenience, we use Ninja as build tool (provided via the ``buildtools`` module).

METIS/ParMETIS:

.. code-block:: bash

    wget https://deb.debian.org/debian/pool/non-free/p/parmetis/parmetis_4.0.3.orig.tar.gz
    tar -xvf parmetis_4.0.3.orig.tar.gz
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
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -DFORTRAN_SUPPORT=off -GNinja -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    ninja install
    cd ../..

For LUA:

.. code-block:: bash

    wget https://www.lua.org/ftp/lua-5.4.6.tar.gz
    tar -xf lua-5.4.6.tar.gz
    cd lua-5.4.6
    make all install INSTALL_TOP=$SEISSOL_PREFIX MYCFLAGS="-fPIC"
    cd ..

For easi (depending on the former two):

.. code-block:: bash

    git clone --recursive --depth 1 https://github.com/seissol/easi
    mkdir -p easi/build
    cd easi/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja -DASAGI=ON -DLUA=ON -DIMPALAJIT=OFF -DEASICUBE=OFF -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())") -DPYTHON_BINDINGS=ON
    ninja install
    cd ../..

For libxsmm (note that we need 1.17 sharp; the latest main will not work as intended with the generator):

.. code-block:: bash

    git clone --branch 1.17 --depth 1 https://github.com/libxsmm/libxsmm
    cd libxsmm
    make generator
    cp bin/libxsmm_gemm_generator $SEISSOL_PREFIX/bin
    cd ..

In case there are problems with using libxsmm, you can also consider using only PSpaMM instead; at a tiny performance penalty.

Compiling SeisSol
~~~~~~~~~~~~~~~~~

Finally, it's time to clone SeisSol and build it.

However, we need to apply a small hotfix here, since the Cray compiler environment does not work with AdaptiveCpp (it causes problems with finding MPI, the filesystem headers etc.). As a workaround, we compile SeisSol with ``amdclang`` directly, and add the necessary flags from the Cray environment as compiler flags (that can be done by ``CC --cray-print-opts=all``, the same with ``cc`` and ``ftn``).

In total, we get the following:

.. code-block:: bash

    git clone --recursive https://github.com/SeisSol/SeisSol.git seissol
    mkdir -p seissol/build
    cd seissol/build
    CC=amdclang CXX=amdclang++ CFLAGS=$(cc --cray-print-opts=all) CXXFLAGS=$(CC --cray-print-opts=all) cmake .. -DPython3_EXECUTABLE=$(which python3.11) -GNinja -DPRECISION=single -DDEVICE_BACKEND=hip -DDEVICE_ARCH=gfx90a -DHOST_ARCH=milan -DORDER=5 -DASAGI=ON -DNUMA_AWARE_PINNING=ON -DUSE_GRAPH_CAPTURING=ON -DDR_QUAD_RULE=dunavant -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
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

    # see https://github.com/SeisSol/SeisSol/issues/1458 for why the --unbuffered
    srun --unbuffered --cpu-bind=mask_cpu:${CPU_BIND} seissol-launch SeisSol_Release_sgfx90a_hip_6_elastic parameters.par
