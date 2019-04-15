Compilation
===========

In order to run SeisSol, you need to first install:

-  Python (use version 2.7 and above)
-  SCons (for instructions see below)
-  hdf5 (for instructions see below)
-  netcdf (C-Release) (for instructions see below)
-  Intel compiler (icc, icpc, ifort) or GCC (gcc, g++, gfortran)
-  Some MPI implementation
-  ParMETIS for partitioning
-  libxsmm (libxsmm\_gemm\_generator) for small matrix multiplications
-  CMake for compiling ImpalaJIT

Inital Adjustments to .bashrc
-----------------------------

Create a new text file in your $HOME directory

.. code-block:: bash

  cd $HOME
  touch seissol_env.sh

Open the file with your favorite text editor (vi ./seissol_env.sh) and add following lines:

.. code-block:: bash

  # User specific aliases and functions
  export S3_HOME=$PWD/s3_software
  
  export PATH=$S3_HOME/bin:$PATH
  export LD_LIBRARY_PATH=$S3_HOME/lib:$S3_HOME/lib64:$LD_LIBRARY_PATH  # lookup path for executables
  export LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH  # lookup path for compilers
  export EDITOR=vi
  export SCONS_LIB_DIR=$S3_HOME/lib64/scons-2.2.0/
  export C_INCLUDE_PATH=$S3_HOME/include:$C_INCLUDE_PATH
  export CPLUS_INCLUDE_PATH=$S3_HOME/include:$CPLUS_INCLUDE_PATH
  export PKG_CONFIG_PATH=$S3_HOME/lib/pkgconfig

  export FC=mpiifort  # or mpifort
  export CXX=mpiicpc  # or mpicxx
  export CC=mpiicc    # or mpicc

  ######  ParMetis library necessary (Optional) ##############
  export PARMETIS_BASE=$S3_HOME
  export PARMETIS_LIBDIR=$S3_HOME/lib


Add the following to your .bashrc (vi ~/.bashrc).

.. code-block:: bash

  source /opt/intel/compiler/VERSION/bin/compilervars.sh intel64 # if you use intel compilers
  source $HOME/seissol_env.sh


.. code-block:: bash

  source ./seissol_env.sh  # update your enviroment variable 

  mkdir s3_software
  cd ./s3_software


Installing SCons
----------------

.. code-block:: bash

  wget http://prdownloads.sourceforge.net/scons/scons-2.2.0.tar.gz
  tar -xvf ./scons-2.2.0.tar.gz
  cd scons-2.2.0
  python setup.py install --prefix=$S3_HOME
  cd ..


Installing HDF5
---------------

.. code-block:: bash

  git clone https://github.com/mortenpi/hdf5.git
  cd hdf5
  git checkout 8a275ab7831002f3e  # take a stable version (optional)
  ./configure --prefix=$S3_HOME --enable-parallel --with-zlib --disable-shared --enable-fortran
  make -j4
  make install
  cd ..


Installing netCDF
-----------------

.. code-block:: bash

  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz
  tar -xaf ./netcdf-*.tar.gz
  cd netcdf-4.4.1.1
  CPPFLAGS="-I$S3_HOME/include -fPIC" ./configure --prefix=$S3_HOME --enable-shared=no
  make -j4
  make check # to check correctness of installation (optional)
  make install
  cd ..



Installing Libxsmm
------------------

.. code-block:: bash

  git clone https://github.com/hfp/libxsmm
  cd libxsmm
  git checkout b6de187f832a723295a  # take a stable version (optional) 
  make generator
  cp ./bin/libxsmm_gemm_generator $S3_HOME/bin/
  cd ..


Installing Metis (Optional:PUML mesh format)
--------------------------------------------

.. code-block:: bash

  wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
  tar -xvf ./metis-5.1.0.tar.gz
  cd metis-5.1.0
  make config cc=$CC cxx=$CXX prefix=$S3_HOME
  make install
  cd ..


Installing ParMetis (Optional:PUML mesh format)
-----------------------------------------------

.. code-block:: bash

  wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
  tar -xvf ./parmetis-4.0.3.tar.gz
  cd parmetis-4.0.3
  make config cc=$CC cxx=$CXX prefix=$S3_HOME
  make install
  cd ..


Installing ASAGI (Optional)
---------------------------

.. code-block:: bash

  git clone https://github.com/TUM-I5/ASAGI.git
  cd ASAGI

  git clone https://github.com/TUM-I5/utils.git
  mkdir build
  cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=$S3_HOME
  make -j4
  make install
  cd ../..


For details, see section :ref:`Installing ASAGI <installing_ASAGI>`.

.. _compiling-seissol:

Compiling SeisSol
-----------------

Get the latest version of SeisSol on git by cloning the whole repository
including all submodules:

.. code-block:: bash

   git clone https://github.com/SeisSol/SeisSol.git
   cd SeisSol
   git submodule update --init

Add the following build variables to the file
build/options/supermuc_mac_cluster.py

.. code-block:: python

   compileMode='release' 
   parallelization='hybrid' 
   arch='$ARCH' 
   order='$ORDER' 
   generatedKernels = 'yes'
   compiler = 'intel' # or gcc
   logLevel = 'info'

   netcdf='yes' 
   hdf5='yes'
   metis = 'yes' #  additionally for puml mesh format
   asagi = 'yes' #  optional for ASAGI


| with: 
| compileMode - release / relWithDebInfo/ debug
| parallelization - omp/ mpi / hybrid (mpi/openmp)
| logLevel - info/ debug, warning or error 
| ARCH - target architecture 
| ORDER - convergence order (=max polynomial order +1)
| generatedKernels - yes/no

Get your executable with

.. code-block:: bash

   scons -j 4 buildVariablesFile=build/options/supermuc_mac_cluster.py

   # choose a correct executable file in case if you compile multiple
   # versions of Seissol
   ln -s $PWD/build/SeisSol_<version> $S3_HOME/bin/Seissol
   cd $HOME


NOTE: SCons will try to detect the correct MPI wrappers. If this fails,
you can overwrite the detected wrappers with the variables "mpicc",
"mpicxx" and "mpif90".

you can run ``scons -h`` to get some help on options

Please note, this builds the generated kernel version of SeisSols. For
SeisSol classic, please add the generatedKernels=no switch. However this
result in roughly 6X less performance. The classic version won't be
maintained anymore in the near future.

Running SeisSol
---------------

1. Follow the instructions on :ref:`Configuration <Configuration>`.
2. run SeisSol version of interest. To run the example:
   ``SeisSol <PARAMETER_FILE_NAME>.par``

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also :ref:`A first example <a_first_example>`.
