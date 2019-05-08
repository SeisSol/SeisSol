Supermuc-NG
~~~~~~~~~~~~

Prerequisites: initialization GitHub as SuperMuc Phase2 on the last page.

1. git clone git@github:SeisSol/SeisSol.git  

2. Load module. Could add these lines to .bashrc:

::

  ##### module load for SeisSol
  module load scons cmake/3.6 python/2.7_intel slurm_setup
  module load parmetis/4.0.3-intel-impi-i64-r64 metis/5.1.0-intel-i64-r64
  module load netcdf/4.6.1-intel-impi-hdf5v1.8-parallel hdf5/1.8.20-intel-impi
  module switch spack/staging/19.1 spack/master
  module switch devEnv/Intel devEnv/Intel/2019
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


3. Install libxsmm and ASAGI (Optional)

Install libxsmm:

::

  $ git clone git@github:hfp/libxsmm.git
  $ Cd libxsmm
  $ make generator
  (then add path-to-libxsmm/bin to PATH in bashrc)

Install ASAGI:

::

  $ git clone git@github:TUM-I5/ASAGI.git
  
then copy fix_submodules.sh here in submodules/utils and run it to get all submodules cloned.

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
  (Know errors: 1.Numa could not found - turn off Numa in the CmakeList.txt. )


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
  zlibDir='$HOME/myLib/ASAGI/build'

  phase=3
  if phase==1:
     arch ='dsnb'
  elif phase==2:
     arch = 'dhsw'
     #commThread ='yes'
  else:
     arch = 'dskx'
     commThread ='yes'

  plasticity='yes'
  #logLevel                    = 'warning'
  logLevel                    = 'warning'
  logLevel0                   = 'info'


5. Submission file for SeisSo on Phase 2:

::

  #!/bin/bash
  # this job command file is called submit.cmd
  #@ energy_policy_tag = <user id>_etag
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
  #@ notify_user=<your email address>
  #@ queue
  
  . /etc/profile
  . /etc/profile.d/modules.sh

  export PARMETIS_BASE='/lrz/sys/libraries/parmetis/4.0.2/ibmmpi'
  export PARMETIS_LIBDIR='/lrz/sys/libraries/parmetis/4.0.2/ibmmpi/lib'

  export MP_SINGLE_THREAD=yes
  export OMP_NUM_THREADS=16

  export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

  #module load mpi_pinning/hybrid_blocked
  #export LD_LIBRARY_PATH=~/

  # ############## dsnb for phase 1 and dhsw for phase 2 ###########################
  cd <work-directory>
  poe ./SeisSol_release_generatedKernels_dhsw_hybrid_none_9_4  parameters.par
  echo "JOB is run"

