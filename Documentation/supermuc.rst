.. _compile_run_supermuc:

Compiling and running SeisSol on Supermuc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Accessing github from SuperMUC
==============================

SuperMUC restricts access to outside sources and thus does not allow connections to https servers. 
Nevetheless, github can be used if remote port forwarding is correctly set.
Here, we described the procedure to setup such port forwarding.


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
  module load scons cmake/3.6 python/3.6_intel slurm_setup
  module load parmetis/4.0.3-intel-impi-i64-r64 metis/5.1.0-intel-i64-r64
  module load netcdf/4.6.1-intel-impi-hdf5v1.8-parallel hdf5/1.8.20-intel-impi-threadsafe
  module load libszip/2.1.1

  ####### universal setup for SeisSol
  export PATH=~/local/bin:$PATH
  export PKG_CONFIG_PATH=~/local/lib/pkgconfig/:$PKG_CONFIG_PATH
  export LD_LIBRARY_PATH=~/local/lib:$PARMETIS_LIB:$METIS_LIB:$NETCDF_BASE/lib:$HDF5_BASE/lib:$LD_LIBRARY_PATH
  export CPATH=~/local/include:$PARMETIS_INC:$METIS_INC:$NETCDF_BASE/include:$HDF5_BASE/include:$CPATH
  export LIBRARY_PATH=~/local/lib:$PARMETIS_LIB:$METIS_LIB:$NETCDF_BASE/lib:$HDF5_BASE/lib:$LIBRARY_PATH
    
  ####  local setup for SeisSol. I installed libxsmm and ASAGI in my own folder (/dss/dsshome1/02/di52lak2/myLib). 
  export PATH=/dss/dsshome1/02/di52lak2/myLib/libxsmm/bin:$PATH
  export PKG_CONFIG_PATH=/dss/dsshome1/02/di52lak2/myLib/ASAGI/build/lib/pkgconfig:$PKG_CONFIG_PATH
  export LD_LIBRARY_PATH=/dss/dsshome1/02/di52lak2/myLib/ASAGI/build/lib:$LD_LIBRARY_PATH


3. Install libxsmm and ASAGI

| See :ref:`installing_libxsmm` and :ref:`installing_ASAGI`. 
| Note that on project pr63qo, we already installed and shared these library (no need to install).
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

  $ export FC=mpif90
  $ export CXX=mpiCC
  $ export CC=mpicc

  $ make build
  $ Cd build
  $ CMAKE_PREFIX_PATH=$NETCDF_BASE
  $ cmake ../ -DSHARED_LIB=no -DSTATIC_LIB=yes -DNONUMA=on -DCMAKE_INSTALL_PREFIX=$HOME/<folder-to-ASAGI>/build/ 
  $ make
  $ make install
  (Know errors: 1.Numa could not found - turn off Numa by -DNONUMA=on . )


4. Copy the SeisSol configuration to a file e.g. supermuc_ng.py

::

  import os
  # build options
  compileMode                 = 'release'
  #compileMode                 = 'relWithDebInfo'
  #compileMode                 = 'debug'
  parallelization             = 'hybrid'
  #parallelization             = 'mpi'
  generatedKernels            = 'yes'
  #measureNodeLevelPerformance = 'none'
  useExecutionEnvironment     = 'yes'
  order = 4
  equations='elastic'
  #equations = 'viscoelastic2'
  #numberOfMechanisms = 3
  # machine dependent options
  #compiler='gcc'
  cppCompiler          = 'mpiCC'
  fortranCompiler      = 'mpif90'

  netcdf='yes'
  hdf5='yes'
  metis='yes'
  netcdfDir=os.environ['NETCDF_BASE']
  hdf5Dir=os.environ['HDF5_BASE']

  # ASAGI folder need to be verified.
  asagi='yes'
  zlibDir='/dss/dsshome1/02/di52lak2/myLib/ASAGI/build/lib'

  phase=3 # for Supermuc-NG
  if phase==1:
     arch ='dsnb'
  elif phase==2:
     arch = 'dhsw'
     #commThread ='yes'
  else:
     arch = 'dskx'
     commThread ='yes'

  plasticity='no'
  #logLevel                    = 'warning'
  logLevel                    = 'warning'
  logLevel0                   = 'info'


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
  srun ./SeisSol_release_generatedKernels_dskx_hybrid_none_9_4 parameters.par

  


 
Supermuc-2
==========

1. Load modules

You can create a folder ~/.modules and copy these to ~/.modules/bash (Must use intel/17.0)
:: 

  module load python/3.5_intel
  module load scons/3.0.1
  module unload netcdf
  module load netcdf/mpi
  module load hdf5/mpi/1.8.18
  module unload intel
  module load intel/17.0
  module load gcc
  module load cmake
  module load szip
  Module load parmetis/4.0

2. Add these lines to .bashrc (there are shared libs under /home/hpc/pr63qo/di52lak/ but you can install by yourself)
::

  ## need installation before and added here ######
  export PATH = $PATH:$HOME/bin:/home/hpc/pr63qo/di52lak/software/libxsmm-master/bin
  export PKG_CONFIG_PATH =/home/hpc/pr63qo/di52lak/software/ASAGI/build/lib/pkgconfig
  export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:/home/hpc/pr63qo/di52lak/software/ASAGI/build/lib

  ## existing lib in supermuc
  export LD_LIBRARY_PATH = /usr/lib64:/lib64:/lib:$LD_LIBRARY_PATH

  ######  parmetis library necessary  ##############
  export PARMETIS_BASE='/lrz/sys/libraries/parmetis/4.0.2/ibmmpi'
  export PARMETIS_LIBDIR='/lrz/sys/libraries/parmetis/4.0.2/ibmmpi/lib'
  export PATH=$PATH:/lrz/sys/libraries/hdf5/1.8.14/ibmmpi_poe1.4_15.0.5/bin
 

3. Build variable file — updated on July 2018

Copy this to a supermuc_hw.py file in SeisSol/:
::

  import os
  # build options
  compileMode = 'release' # or relWithDebInfo or debug
  generatedKernels = 'yes'
  arch = 'dhsw'  # use 'dsnb' for SuperMUC phase 1 or use 'dhsw' for SuperMUC phase 2
  parallelization = 'hybrid'
  order = '4' # valid values are 'none', '2', '3', '4', '5', '6', '7', and '8'.
  equations = 'elastic' # valid values are 'elastic', 'viscoelastic', 'viscoelastic2'
  plasticity = 'no' # start with elastic at the beginning.

  useExecutionEnvironment = 'yes'
  logLevel = 'warning'
  logLevel0 = 'info'

  netcdf = 'yes'
  hdf5 = 'yes'
  metis = 'yes'
  netcdfDir=os.environ['NETCDF_BASE']
  hdf5Dir=os.environ['HDF5_BASE']
  metisDir = '/lrz/sys/libraries/parmetis/4.0.2/ibmmpi'

  asagi = 'yes’
  zlibDir='/home/hpc/pr63po/di52lak/software/ASAGI/build/lib/'

  # Put a 'yes' here on Phase 2 and a 'no' on Phase 1
  commThread = 'no'
  # If you put a 'yes' for the last option on Phase 2, it is vital that your environment settings are correct, otherwise your performance will be bad.


4. compile SeisSol as:

::

  scons buildVariablesFile=supermuc_hw.py
  
  
5. Submit job on Phase 2. Here is an example:

::

  #!/bin/bash
  # this job command file is called submit.cmd
  #@ energy_policy_tag = <account id>_etag
  #@ minimize_time_to_solution = yes
  #@ wall_clock_limit = 12:00:00

  #@ job_name = <job name>
  #@ class = micro
  #@ island_count=1
  ## #@ input= job.$(schedd_host).$(jobid).in
  #@ output= job.$(schedd_host).$(jobid).out
  #@ error= job.$(schedd_host).$(jobid).err
  #@ job_type= parallel 
  #@ node= 7
  #@ tasks_per_node= 1
  ## #@ total_tasks= 512
  #@ network.MPI = sn_all,not_shared,us
  #@ notification=always
  #@ notify_user=dli@geophysik.uni-muenchen.de
  #@ queue
  . /etc/profile
  . /etc/profile.d/modules.sh

  export PARMETIS_BASE='/lrz/sys/libraries/parmetis/4.0.2/ibmmpi'
  export PARMETIS_LIBDIR='/lrz/sys/libraries/parmetis/4.0.2/ibmmpi/lib'

  export MP_SINGLE_THREAD=yes
  export OMP_NUM_THREADS=28
  export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS


  # ############## dsnb for phase 1 and dhsw for phase 2 ###########################
  cd <working directory>
  poe ./SeisSol_release_generatedKernels_dhsw_hybrid_none_9_4 parameters.par
  echo "JOB is run"

see the enviroment variables section concerning :ref:`environement_variables_supermuc_phase_2` for more details.

