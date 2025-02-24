..
  SPDX-FileCopyrightText: 2021-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_marconi:

Marconi 100
===========

Marconi 100 is a distributed multi-GPU HPC system equipped with 4 Nvidia V100 GPUs
and 2 IBM Power9 CPUs per node. This architecture usually comes with a pre-installed
CUDA-Aware Spectrum-MPI. However, SeisSol cannot operate with Spectrum-MPI because
of GPU memory leaks caused by this MPI implementation. This part of documentation
describes how to setup and configure OpenMPI 4.1.x together with UCX 1.10.x to operate
on IBM Power and similar HPC systems.

Installing Main Libraries and Packages
---------------------------------------

1. Create the following launch file:

::

  $ touch $HOME/launch.sh
  $ cat $HOME/launch.sh
  #!/bin/bash

  VERSION="openmpi4-1-1_ucx1-10-1"
  export SEISSOL_INSTALL=$HOME/usr/local/${VERSION}

  export PATH=$SEISSOL_INSTALL/bin:$PATH
  export LIBRARY_PATH=$SEISSOL_INSTALL/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$SEISSOL_INSTALL/lib:$LD_LIBRARY_PATH
  export PKG_CONFIG_PATH=$SEISSOL_INSTALL/lib/pkgconfig:$PKG_CONFIG_PATH
  export CMAKE_PREFIX_PATH=$SEISSOL_INSTALL
  export CPATH=$SEISSOL_INSTALL/include:$CPATH
  export CPLUS_INCLUDE_PATH=$SEISSOL_INSTALL/include:$CPLUS_INCLUDE_PATH
  export C_INCLUDE_PATH=$SEISSOL_INSTALL/include:$C_INCLUDE_PATH

2. Load basic modules, source launch file and create an install directory:

::

  $ source ./launch.sh
  $ mkdir -p $SEISSOL_INSTALL

  $ module load python/3.8.2
  $ module load profile/candidate
  $ module load gnu/11.2.0
  $ module load cmake/3.20.0
  $ module load cuda/11.6

3. Create any directory where you are going to configure and build libraries and packages. For example,

::

  mkdir -p $HOME/Download
  cd $HOME/Download

4. Install SYCL:

See section :ref:`Installing SYCL <installing_SYCL>`. Note, you will need to adjust
install-prefixes for both `LLVM`, `Boost` and `hipSYCL` to point to the content of `SEISSOL_INSTALL` environment
variable.

5. Install hwloc:

::

  $ wget http://www.open-mpi.org/software/hwloc/v2.7/downloads/hwloc-2.7.1.tar.gz && \
  $ tar -xvf ./hwloc-2.7.1.tar.gz
  $ cd hwloc-2.7.1 && \
  $ CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran) ./configure --prefix=${SEISSOL_INSTALL} \
    --disable-opencl --disable-cairo \
    --disable-nvml --disable-gl \
    --enable-cuda \
    --with-cuda=/cineca/prod/opt/compilers/cuda/11.6/none/bin/nvcc \
    --disable-libudev --enable-shared

  $ make -j
  $ make install
  $ cd ..

6. Install UCX:

Note, CINECAs' system administrators manually added `CUDA_CFLAGS` env. variable to all `cuda` modules. This prevents
compilation of many versions of UCX and OpenMPI. Please, disable it during SeisSol's installation.

::

  $ wget https://github.com/openucx/ucx/releases/download/v1.10.0/ucx-1.10.0.tar.gz
  $ tar -xvf ucx-1.10.0.tar.gz
  $ unset CUDA_CFLAGS
  $ mkdir -p ucx-1.10.0/build && cd ucx-1.10.0/build

  $ CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran) \
  ../configure \
  --prefix=$SEISSOL_INSTALL \
  --build=powerpc64le-redhat-linux-gnu \
  --host=powerpc64le-redhat-linux-gnu \
  --disable-logging --disable-debug \
  --disable-assertions --disable-params-check \
  --disable-params-check --enable-optimizations \
  --disable-assertions --disable-logging --with-pic \
  --without-java \
  --enable-mt \
  --with-cuda=/cineca/prod/opt/compilers/cuda/11.6/none/bin/nvcc \
  --with-gdrcopy \
  --with-knem=/opt/knem-1.1.3.90mlnx1 \
  --without-xpmem

  $ make -j
  $ make install
  $ cd ../..

7. Install OpenMPI:

::

  $ wget https://github.com/open-mpi/ompi/archive/refs/tags/v4.1.4.tar.gz
  $ tar -xvf ./v4.1.4.tar.gz
  $ cd ./ompi-4.1.4
  $ ./autogen.pl
  $ mkdir -p ./build && cd ./build

  $ CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran) \
  CFLAGS="-I/opt/pmix/3.1.5/include" CPPFLAGS="-I/opt/pmix/3.1.5/include" \
  ../configure \
  --prefix=$SEISSOL_INSTALL \
  --with-memory-manager=none \
  --enable-static=yes \
  --enable-shared \
  --with-slurm \
  --with-pmix=/opt/pmix/3.1.5 \
  --with-ucx=$SEISSOL_INSTALL \
  --with-libevent=/usr \
  --with-hwloc=${SEISSOL_INSTALL} \
  --with-verbs \
  --enable-mpirun-prefix-by-default \
  --with-platform=/cineca/prod/build/compilers/openmpi/4.0.3/gnu--8.4.0/BA_WORK/openmpi-4.0.3/contrib/platform/mellanox/optimized \
  --with-hcoll=/opt/mellanox/hcoll \
  --with-cuda=/cineca/prod/opt/compilers/cuda/11.6/none/bin/nvcc \
  --with-knem=/opt/knem-1.1.3.90mlnx1 \
  --without-xpmem

  $ make -j
  $ make install
  $ cd ../..

8. Install HDF5:

::

  $ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz
  $ tar -xvf ./hdf5-1.10.5.tar.gz
  $ cd hdf5-1.10.5
  $ ./autogen.sh
  $ mkdir build && cd build

  $ CFLAGS="-fPIC ${CFLAGS}" CC=mpicc CXX=mpicxx FC=mpif90 \
  ../configure \
  --prefix=$SEISSOL_INSTALL \
  --build=powerpc64le-redhat-linux-gnu \
  --host=powerpc64le-redhat-linux-gnu \
  --enable-parallel --with-zlib --disable-shared \
  --enable-fortran

  $ make -j
  $ make install
  $ cd ../..

9. Installing netCDF:

::

  $ wget https://syncandshare.lrz.de/dl/fiJNAokgbe2vNU66Ru17DAjT/netcdf-4.6.1.tar.gz
  $ tar -xvf ./netcdf-4.6.1.tar.gz
  $ cd hdf5-1.10.5
  $ ./autogen.sh

  $ CFLAGS="-fPIC ${CFLAGS}" CC=h5pcc \
  ./configure \
  --prefix=$SEISSOL_INSTALL \
  --build=powerpc64le-redhat-linux-gnu \
  --host=powerpc64le-redhat-linux-gnu \
  --enable-shared=no \
  --disable-dap

  $ make -j
  $ make install
  $ cd ..

10. Installing ParMetis:

::

  $ https://ftp.mcs.anl.gov/pub/pdetools/spack-pkgs/parmetis-4.0.3.tar.gz
  $ tar -xvf ./parmetis-4.0.3.tar.gz
  $ cd parmetis-4.0.3
  #edit ./metis/include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
  $ make config cc=mpicc cxx=mpiCC prefix=$SEISSOL_INSTALL
  $ make install
  $ cp build/Linux-ppc64le/libmetis/libmetis.a $SEISSOL_INSTALL/lib
  $ cp metis/include/metis.h $SEISSOL_INSTALL/include
  $ cd ..

11. Install GemmForge and ChainForge. Please, follow steps described :ref:`here <gemmforge_installation>`.

12. Install easi with LUA backend:

::

  # yaml-cpp
  $ wget https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.7.0.tar.gz
  $ tar -xvf yaml-cpp-0.7.0.tar.gz
  $ cd yaml-cpp-yaml-cpp-0.7.0
  $ sed -i 's/$<${not-msvc}/#$<${not-msvc}/g' ./CMakeLists.txt
  $ mkdir -p build && cd build
  $ cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_INSTALL \
    -DYAML_BUILD_SHARED_LIBS=ON \
    -DBUILD_TESTING=OFF \
    -DBUILD_MOCK=OFF \
    -DBUILD_GMOCK=OFF
  $ make -j4 && make install
  $ cd ../..

  # LUA
  $ wget https://www.lua.org/ftp/lua-5.3.6.tar.gz
  $ tar -xzvf lua-5.3.6.tar.gz
  $ cd lua-5.3.6
  $ make linux CC=mpicc
  $ make local
  $ cp -r install/* $SEISSOL_INSTALL
  $ cd ..

  # easi
  $ git clone https://github.com/SeisSol/easi.git
  $ cd easi
  $ mkdir -p build && cd build
  $ CC=mpicc CXX=mpicxx FC=mpifort cmake .. \
    -DASAGI=OFF \
    -DLUA=ON \
    -DIMPALAJIT=OF \
    -DCMAKE_INSTALL_PREFIX=$SEISSOL_INSTALL
  $ make -j4 && make install
  $ cd ../..

13. Install Eigen3:

::

  $ git clone https://gitlab.com/libeigen/eigen.git
  $ mkdir -p eigen/build && cd eigen/build
  $ CXX=g++ CC=gcc FC=gfortran cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_INSTALL
  $ make -j
  $ make install
  $ cd ../..

14. Install SeisSol:

::

  $ module load cuda/11.6
  $ git clone --recurse-submodules https://github.com/SeisSol/SeisSol.git
  $ cd SeisSol
  $ mkdir build && cd build

  $ CC=mpicc CXX=mpicxx FC=mpifort cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DDEVICE_BACKEND=cuda \
  -DDEVICE_ARCH=sm_70 \
  -DHOST_ARCH=power9 \
  -DPRECISION=single

  $ make -j

15. Run SeisSol-proxy as a sanity check:

::

  ./launch ./SeisSol_proxy_Release_ssm70_cuda_6_elastic 100000 100 all


Running SeisSol
---------------

As discussed :ref:`here <gpu_process_pinning>`, process pinning is important for SeisSol GPU version.
IBM Power9 is an example of RISC architecture designed with with 4-way hyperthreading and 8 cores per CPU.
In total, each node of Marconi 100 can run 256 threads. By and large process pinning needs a special
care on such architectures because some libraries have different meanings of cores and threads.

Below you can see an example of a *batch script* with parameters resulting in an optimal process pinning.
Note that each node of Marconi 100 has 2 Mellanox network cards i.e., each per NUMA domain. In this example,
we enforce UCX to utilize both. Moreover, we reserve one 1 core for each MPI process for SeisSol communication thread.

In this particular case it is not necessary to provide a number of processes after **mpirun** because OpenMPI was compiled
with PMIX (see step 5).

Please, do not forget to launch SeisSol via `launch` bash script generated with CMake during SeisSol's configuration.

::

  $ cat ./job.sh
  #!/bin/bash
  #SBATCH --account=<you account>
  #SBATCH --partition=m100_usr_prod
  #SBATCH --qos=m100_qos_dbg
  #SBATCH --time <time>
  #SBATCH --nodes <number of nodes>
  #SBATCH --ntasks-per-node=4
  #SBATCH --cpus-per-task=32
  #SBATCH --gres=gpu:4
  #SBATCH --gpu-bind=closest
  #SBATCH --mem=161070
  #SBATCH --job-name=<your job name>
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=<user_email>
  #SBATCH --output=seissol.out
  #SBATCH --error=seissol.err
  #SBATCH --exclusive

  NUM_CORES=$(expr $SLURM_CPUS_PER_TASK / 4)
  NUM_COMPUTE_CORES=$(expr $NUM_CORES - 1)

  export OMP_NUM_THREADS=$NUM_COMPUTE_CORES
  export OMP_PLACES="cores($NUM_COMPUTE_CORES)"
  export PROC_BIND=spread

  export DEVICE_STACK_MEM_SIZE=1.5
  export UCX_MEMTYPE_CACHE=n

  mpirun --report-bindings --map-by ppr:$SLURM_NTASKS_PER_NODE:node:pe=$NUM_CORES \
  -x UCX_MAX_EAGER_RAILS=2 -x UCX_MAX_RNDV_RAILS=2 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_1:1 \
  -x UCX_MEM_MMAP_HOOK_MODE=none \
  ./launch ./SeisSol_Release_ssm70_cuda_6_elastic ./parameters.par



Troubleshooting OpenMPI and UCX
-------------------------------

1. OpenMPI and UCX provide users with utilities which show how these packages were configured and installed.
It is **ompi_info** for former and **ucx_info -b** for latter.

2. It may be required to switch off UCX memory caching because it can lead to run-time failures of UCX.
One can achieve this by setting the following environment variable:

::

  $ export UCX_MEMTYPE_CACHE=n

3. One can enable a launch time information from OpenMPI and UCX by setting the following parameters after **mpirun**.

::

  --mca pml_base_verbose 10 --mca mtl_base_verbose 10 -x OMPI_MCA_pml_ucx_verbose=10

4. We recommend to login into a compute node and execute **ucx_info -d**  command if you need to get information
about all available network devices. This will help you to retrieve exact names of network devices e.g., *mlx5_0:1, mlx5_1:1, etc*.
