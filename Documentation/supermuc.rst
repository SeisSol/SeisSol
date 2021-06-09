.. _compile_run_supermuc:


Compiling and running SeisSol on Supercomputers
===============================================

Supermuc
--------

Accessing github from SuperMUC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SuperMUC restricts access to outside sources and thus does not allow connections to https servers. 
Nevertheless, GitHub can be used if remote port forwarding is correctly set.
Here, we described the procedure to set up such port forwarding.


1. On your local machine, add to your ~/.ssh/config the following lines:

::

  Host supermuNG
     Hostname skx.supermuc.lrz.de
     User <Your Login>    
     RemoteForward ddddd github.com:22

where ddddd is an arbitrary 5-digital port number.

2. ssh supermucNG to login to supermuc. Then add the following lines to the ~/.ssh/config:

:: 

  Host github.com
     HostName localhost
     User git
     Port ddddd
    
With ddddd the same port number as before.

4. Create SSH key by typing (use a non-empty passphrase, not too long as you will need to type it often)

::

  ssh-keygen -t rsa 

5. Go to https://github.com/settings/ssh, add a new SSH key, pasting the public key you just created on supermuc  ~/.ssh/id_rsa.pub. 
Logout of supermuc and log back in (ssh supermucNG). You should now be able to clone SeisSol including the submodules using:


::

  git clone git@github.com:SeisSol/SeisSol.git
  cd SeisSol
  git submodule update --init

Pay attention to the git clone address ('https://github.com/' replaced by 'git@github.com:'). 
If it works, you will see several lines of ‘cloning ….’.


Supermuc-NG
~~~~~~~~~~~

1. clone SeisSol including the submodules using 

::

  git clone git@github.com:SeisSol/SeisSol.git
  cd SeisSol
  git submodule update --init
 

2. Load module. Could add these lines to .bashrc (changing the order and adding additionnal modules may prevent a successful compilation):

::

  ##### module load for SeisSol
  module load scons gcc/9 cmake/3.14.4 python/3.6_intel
  module load libszip/2.1.1
  module load parmetis/4.0.3-intel-impi-i64-r64 metis/5.1.0-intel-i64-r64
  module load hdf5/1.8.21-impi-cxx-frt-threadsafe 
  module load netcdf/4.6.1-intel-impi-hdf5v1.8-parallel
  module load numactl

  ####### for pspamm.py
  export PATH=~/bin:$PATH
  
  ####  local setup for SeisSol. 
  export PATH=/hppfs/work/pr63qo/di73yeq4/myLibs/libxsmm/bin:$PATH
  export PKG_CONFIG_PATH=/hppfs/work/pr63qo/di73yeq4/myLibs/ASAGI/build/lib/pkgconfig:$PKG_CONFIG_PATH
  export LD_LIBRARY_PATH=/hppfs/work/pr63qo/di73yeq4/myLibs/ASAGI/build/lib:$LD_LIBRARY_PATH


3. Install libxsmm, PSpaMM and ASAGI

See :ref:`installing_libxsmm`, :ref:`installing_pspamm` and :ref:`installing_ASAGI`. 
Note that on project pr63qo, we already installed and shared libxsmm and ASAGI (but not pspamm).
The compiled libs are in /hppfs/work/pr63qo/di73yeq4/myLibs/xxxx/build with xxxx=ASAGI or libxsmm.
If you need to compile ASAGI, first clone ASAGI with:

.. code-block:: bash

  git clone git@github.com:TUM-I5/ASAGI
  cd ASAGI
  git submodule update --init
 
set compiler options, run cmake, and compile with:

::

  export FC=mpif90
  export CXX=mpiCC
  export CC=mpicc

  mkdir build && cd build
  CMAKE_PREFIX_PATH=$NETCDF_BASE
  cmake ../ -DSHARED_LIB=no -DSTATIC_LIB=yes -DNONUMA=on -DCMAKE_INSTALL_PREFIX=$HOME/<folder-to-ASAGI>/build/ 
  make -j 48
  make install
  (Know errors: 1.Numa could not found - turn off Numa by adding -DNONUMA=on . )


4. Install SeisSol with cmake, e.g. with (more options with ccmake)

::

   mkdir build-release && cd build-release
   CC=mpicc CXX=mpiCC FC=mpif90  cmake -DCOMMTHREAD=ON -DNUMA_AWARE_PINNING=ON -DASAGI=ON -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=skx -DPRECISION=single -DORDER=4 -DCMAKE_INSTALL_PREFIX=$(pwd)/build-release -DGEMM_TOOLS_LIST=LIBXSMM,PSpaMM -DPSpaMM_PROGRAM=~/bin/pspamm.py ..
   make -j 48

5. Submission file for SeisSol on NG:

::

  #!/bin/bash
  # Job Name and Files (also --job-name)

  #SBATCH -J <job name>
  #Output and error (also --output, --error):
  #SBATCH -o ./%j.%x.out
  #SBATCH -e ./%j.%x.err

  #Initial working directory (also --chdir):
  #SBATCH --workdir=<work directory>

  #Notification and type
  #SBATCH --mail-type=END
  #SBATCH --mail-user=<your email address>

  # Wall clock limit:
  #SBATCH --time=03:00:00
  #SBATCH --no-requeue

  #Setup of execution environment
  #SBATCH --export=ALL
  #SBATCH --account=<project id>
  #constraints are optional
  #--constraint="scratch&work"
  #SBATCH --partition=general

  #Number of nodes and MPI tasks per node:
  #SBATCH --nodes=40
  #SBATCH --ntasks-per-node=1
  module load slurm_setup
  
  #Run the program:
  export MP_SINGLE_THREAD=no
  unset KMP_AFFINITY
  export OMP_NUM_THREADS=94
  export OMP_PLACES="cores(47)"

  export XDMFWRITER_ALIGNMENT=8388608
  export XDMFWRITER_BLOCK_SIZE=8388608
  export SC_CHECKPOINT_ALIGNMENT=8388608

  export SEISSOL_CHECKPOINT_ALIGNMENT=8388608
  export SEISSOL_CHECKPOINT_DIRECT=1
  export ASYNC_MODE=THREAD
  export ASYNC_BUFFER_ALIGNMENT=8388608
  source /etc/profile.d/modules.sh

  echo $SLURM_NTASKS
  ulimit -Ss 2097152
  mpiexec -n $SLURM_NTASKS SeisSol_Release_sskx_4_elastic parameters.par


Marconi 100
-----------

Marconi 100 is a distributed multi-GPU HPC system equipped with 4 Nvidia V100 GPUs
and 2 IBM Power9 CPUs per node. This architecture usually comes with a pre-installed
CUDA-Aware Spectrum-MPI. However, SeisSol cannot operate with Spectrum-MPI because 
of GPU memory leaks caused by this MPI implementation. This part of documentation
describes how to setup and configure OpenMPI 4.1.x together with UCX 1.10.x to operate
on IBM Power and similar HPC systems. 

Installing Main Libraries and Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create the following launch file:

::

  $ touch $HOME/launch.sh
  $ cat $HOME/launch.sh
  #!/bin/bash

  VERSION="openmpi4-1-x_ucx1-10-x"
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

  $ module load gnu/8.4.0
  $ module load cmake/3.20.0
  $ module load python/3.8.2

3. Create any directory where you are going to configure and build libraries and packages. For example, 

::

  mkdir -p $HOME/Download
  cd $HOME/Download

4. Install UCX:

::

  $ git clone https://github.com/openucx/ucx.git
  $ cd ucx
  $ git checkout v1.10.x
  $ ./autogen.sh
  $ mkdir build && cd build

  $ CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran) \
  ../contrib/configure-opt \
  --prefix=$SEISSOL_INSTALL \
  --build=powerpc64le-redhat-linux-gnu \
  --host=powerpc64le-redhat-linux-gnu \
  --disable-logging --disable-debug \
  --disable-assertions --disable-params-check \
  --disable-params-check --enable-optimizations \
  --disable-assertions --disable-logging --with-pic \
  --without-java \
  --enable-mt \
  --with-cuda=/cineca/prod/opt/compilers/cuda/11.0/none \
  --with-gdrcopy \
  --with-knem=/opt/knem-1.1.3.90mlnx1 \
  --without-xpmem

  $ make -j
  $ make install
  $ cd ../..

5. Install OpenMPI:

::

  $ git clone --recursive -b v4.1.x https://github.com/open-mpi/ompi.git
  $ cd ompi
  $ ./autogen.pl
  $ mkdir build && cd build

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
  --with-hwloc=/usr \
  --with-verbs \
  --enable-mpirun-prefix-by-default \
  --with-platform=/cineca/prod/build/compilers/openmpi/4.0.3/gnu--8.4.0/BA_WORK/openmpi-4.0.3/contrib/platform/mellanox/optimized \
  --with-hcoll=/opt/mellanox/hcoll \
  --with-cuda=/cineca/prod/opt/compilers/cuda/11.0/none \
  --with-knem=/opt/knem-1.1.3.90mlnx1 \
  --without-xpmem

  $ make -j
  $ make install
  $ cd ../..

6. Install HDF5:

::

  $ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz
  $ tar -xvf ./hdf5-1.10.5.tar.gz
  $ cd hdf5-1.10.5
  $ ./autogen.sh
  $ mkdir build && cd build

  $ CPPFLAGS="-fPIC ${CPPFLAGS}" CC=mpicc CXX=mpicxx FC=mpif90 \
  ../configure \
  --prefix=$SEISSOL_INSTALL \
  --build=powerpc64le-redhat-linux-gnu \
  --host=powerpc64le-redhat-linux-gnu \
  --enable-parallel --with-zlib --disable-shared \
  --enable-fortran

  $ make -j
  $ make install
  $ cd ../..

7. Installing netCDF:

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

8. Installing ParMetis:

::

  $ wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
  $ tar -xvf ./parmetis-4.0.3.tar.gz
  $ cd parmetis-4.0.3
  #edit ./metis/include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
  $ make config cc=mpicc cxx=mpiCC prefix=$SEISSOL_INSTALL
  $ make install
  $ cp build/Linux-ppc64le/libmetis/libmetis.a $SEISSOL_INSTALL/lib
  $ cp metis/include/metis.h $SEISSOL_INSTALL/include
  $ cd ..

9. Install GemmForge. Please, follow steps described :ref:`here <gemmforge_installation>`. 

10. Install SeisSol:

::

  $ module load cuda/11.0
  $ git clone --recurse-submodules https://github.com/SeisSol/SeisSol.git
  $ cd SeisSol
  $ mkdir build && cd build

  $ CC=mpicc CXX=mpicxx FC=mpifort cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DDEVICE_ARCH=nvidia \
  -DDEVICE_SUB_ARCH=sm_70 \
  -DHOST_ARCH=power9 \
  -DPRECISION=single \
  -DCOMMTHREAD=ON

  $ make -j

11. Run SeisSol-proxy as a sanity check:

::

  ./SeisSol_proxy_Release_snvidia_6_elastic 100000 100 all


Launching SeisSol
~~~~~~~~~~~~~~~~~

As discussed :ref:`here <gpu_process_pinning>`, process pinning is important for SeisSol GPU version.
IBM Power9 is an example of RISC architecture designed with with 4-way hyperthreading and 8 cores per CPU.
In total, each node of Marconi 100 can run 256 threads. By and large process pinning needs a special 
care on such architectures because some libraries have different meanings of cores and threads.

Below you can see an example of a *batch script* with parameters resulting in an optimal process pinning.
Note that each node of Marconi 100 has 2 Mellanox network cards i.e., each per NUMA domain. In this example,
we enforce UCX to utilize both. Moreover, we reserve one 1 core for each MPI process for SeisSol communication thread.

In this particular case it is not necessary to provide a number of processes after **mpirun** because OpenMPI was compiled 
with PMIX (see step 5).
 

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
  #STABCH --gpu-bind=closest
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
  ./SeisSol_Release_snvidia_6_elastic ./parameters.par



Troubleshooting OpenMPI and UCX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
