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

  

