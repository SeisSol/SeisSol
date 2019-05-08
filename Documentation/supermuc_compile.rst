Supermuc Phase2 (hasewell) compilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub Initialization 
======================

1. pick an arbitrary 5-digital number for the port (e.g. 22222 but different for you).

2. On the local computer, go to: > vi ~/.ssh/config (make a .ssh folder first if there is none), copy these lines:

::

  # X11-Forward for all host in our subnet
  # HOST *.geophysik.uni-muenchen.de
  #   ForwardX11 yes 

  HOST *.supermuc.lrz.de
     RemoteForward 22222 github.com:22

  Host supermuc2
     Hostname hw.supermuc.lrz.de
     User <Your Login>    
     RemoteForward 22222 github.com:22

3. ssh supermuc2 to login SuperMuc2. On supermuc go to > vi  ~/.ssh/config and copy this:

:: 

  Host github
     HostName localhost
     User git
     Port 22222
    
(note that port number is the same everywhere!)

4. Create SSH key by typing > ssh-keygen -t rsa and copy all the content in  ~/.ssh/id_rsa.pub  
(make sure the position is ~/.ssh/id_rsa.pub and you can hit return twice to skip the passphrase)
  
5. Go to https://github.com/settings/ssh. Choose to create a new SSH key and paste what you copied here ( The title can be anything you like).
Logout of supermuc and log back in (ssh supermuc2)

::

  git clone git@github:SeisSol/SeisSol.git  


You will see several lines of ‘cloning ….’ means that it works. 
(pay attention to the git clone address. use clone with ssh instead of https. This is the same for anyl other repository such as ASAGI and libxsmm. Note that it should be github instead of github.com)

Now for submodule updates, copy the following to fix_submodules.sh. Put fix_submodules.sh in SeisSol/ and run it.
It clones submodules one by one. 

You can check whether it workd or not by: 

::

  ls submodules/yaml-cpp/*
  

If it is not empty then it worked. This is very useful as every repository need it!

::

  #!/bin/bash                                                                                                            
  git submodule init
  for module_full in `git config -l | awk -F"=" {'print $1'} | grep submodule ` ;
  do
   echo $module_full
   module=$(git config ${module_full} | sed "s/.*github\.com\///")
   if [ "$module" != "true" ]; then
      echo $module
      echo "git config ${module_full} git@github:${module}"
      git config ${module_full} git@github:${module}
   fi
  done
  git submodule update
  
  
Compile SeisSol
================

1. Load modules in Supermuc phase2
You can create a folder ~/.modules and copy these to ~/.modules/bash (Must use intel/17.0)
:: 

  module load python/2.7_anaconda
  module load scons
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
  netcdfDir = '/lrz/sys/libraries/netcdf/4.3.3/intel/ibmmpi_poe1.4_1505'
  hdf5 = 'yes'
  hdf5Dir = '/lrz/sys/libraries/hdf5/1.8.14/ibmmpi_poe1.4_15.0.5'
  metis = 'yes'
  metisDir = '/lrz/sys/libraries/parmetis/4.0.2/ibmmpi'

  asagi = 'yes’
  zlibDir='/home/hpc/pr63po/di52lak/software/ASAGI/build/lib/'

  # Put a 'yes' here on Phase 2 and a 'no' on Phase 1
  commThread = 'yes'
  # If you put a 'yes' for the last option on Phase 2, it is vital that your environment settings are correct, otherwise your performance will be bad.


4. compile SeisSol as:

::

  scons buildVariablesFile=supermuc_hw.py
  
  
5. Submit job on NG. Here is an example for slurm submission file:

::

  #!/bin/bash
  # Job Name and Files (also --job-name)

  #SBATCH -J <job name>
  #Output and error (also --output, --error):
  #SBATCH -o ./%j.%x.out
  #SBATCH -e ./%j.%x.out

  #Initial working directory (also --chdir):
  #SBATCH --workdir=<work directory>

  #Notification and type
  #SBATCH --mail-type=END
  #SBATCH --mail-user=<your email address>

  # Wall clock limit:
  #SBATCH --time=03:00:00
  #SBATCH --no-requeue

  #Setup of execution environment
  #SBATCH --export=NONE
  #SBATCH --account=<project id>
  #constraints are optional
  #--constraint="scratch&work"
  #SBATCH --partition=general

  #Number of nodes and MPI tasks per node:
  #max33 so far, else error
  #SBATCH --nodes=20
  #SBATCH --ntasks-per-node=1
  #SBATCH --cpus-per-task=96
  #Needs specific MPI
  #module switch mpi.intel mpi.intel/2019
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
  export SEISSOL_ASAGI_MPI_MODE=OFF
  source /etc/profile.d/modules.sh

  echo $SLURM_NTASKS
  srun --export=ALL ./SeisSol_release_generatedKernels_dskx_hybrid_none_9_4 parameters.par



