.. _compile_run_supermuc:

Compiling and running SeisSol on Supermuc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Accessing github from SuperMUC
==============================

SuperMUC restricts access to outside sources and thus does not allow connections to https servers. 
Nevertheless, GitHub can be used if remote port forwarding is correctly set.
Here, we described the procedure to set up such port forwarding.


1. Add to you ~/.ssh/config the following lines:

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

4. Create SSH key by typing 

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
===========

1. clone SeisSol including the submodules using 

::

  git clone git@github.com:SeisSol/SeisSol.git
  cd SeisSol
  git submodule update --init
 

2. Load module. Could add these lines to .bashrc:

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

| See :ref:`installing_libxsmm`, :ref:`installing_pspamm` and :ref:`installing_ASAGI`. 
| Note that on project pr63qo, we already installed and shared these libraries (no need to install).
| The compiled libs are in /hppfs/work/pr63qo/di73yeq4/myLibs/xxxx/build with xxxx=ASAGI or libxsmm.
| If you need to compile ASAGI, copy the following to fix_submodules.sh and run it within ASAGI to get submodules/utils cloned.

.. code-block:: bash

  #!/bin/bash                                                                                                            
  git submodule init
  for module_full in `git config -l | awk -F"=" {'print $1'} | grep submodule ` ;
  do
   echo $module_full
   module=$(git config ${module_full} | sed "s/.*github\.com\///")
   if [ "$module" != "true" ]; then
      echo $module
      echo "git config ${module_full} git@github.com:${module}"
      git config ${module_full} git@github.com:${module}
   fi
  done
  git submodule update

set compiler options:

::

  $ export FC=mpiifort
  $ export CXX=mpiicpc
  $ export CC=mpiicc

  $ make build
  $ Cd build
  $ CMAKE_PREFIX_PATH=$NETCDF_BASE
  $ cmake ../ -DSHARED_LIB=no -DSTATIC_LIB=yes -DNONUMA=on -DCMAKE_INSTALL_PREFIX=$HOME/<folder-to-ASAGI>/build/ 
  $ make
  $ make install
  (Know errors: 1.Numa could not found - turn off Numa by adding -DNONUMA=on . )


4. Install SeisSol with cmake, e.g. with (more options with ccmake)

::

   mkdir build-release && cd build-release
   CC=mpiicc CXX=mpiicpc FC=mpiifort  cmake -DCOMMTHREAD=ON -DNUMA_AWARE_PINNING=ON -DASAGI=ON -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=skx -DPRECISION=single -DORDER=4 -DCMAKE_INSTALL_PREFIX=$(pwd)/build-release -DGEMM_TOOLS_LIST=LIBXSMM,PSpaMM -DPSpaMM_PROGRAM=~/bin/pspamm.py ..
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
  srun SeisSol_Release_sskx_4_elastic parameters.par

