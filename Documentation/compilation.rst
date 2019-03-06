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

Add the following lines to your .bashrc (vi ~/.bashrc).

.. code-block:: bash

  source /opt/intel/compiler/VERSION/bin/compilervars.sh intel64

  # User specific aliases and functions
  export PATH=$HOME/bin:$PATH
  export LD_LIBRARY_PATH=$HOME/lib:$HOME/lib64:$LD_LIBRARY_PATH
  export EDITOR=vi
  export SCONS_LIB_DIR=$HOME/lib64/scons-2.2.0/
  export PATH=$HOME/../libxsmm/bin:$PATH
  export C_INCLUDE_PATH=$HOME/include:$C_INCLUDE_PATH 
  export CPLUS_INCLUDE_PATH=$HOME/include:$CPLUS_INCLUDE_PATH

  ######  ParMetis library necessary (Optional) ##############
  export PARMETIS_BASE='path_to_parmetis'
  export PARMETIS_LIBDIR='path_to_parmetis/lib'

Installing SCons
----------------

.. code-block:: bash

  wget http://prdownloads.sourceforge.net/scons/scons-2.2.0.tar.gz
  tar -xaf scons-2.2.0.tar.gz
  cd scons-2.2.0
  python setup.py install --prefix=$HOME
  cd ..


Installing HDF5
---------------

.. code-block:: bash

  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.11/src/hdf5-1.8.11.tar.bz2
  tar -xaf hdf5-1.8.11.tar.bz2
  cd hdf5-1.8.11
  CC=mpiicc FC=mpiifort ./configure --enable-parallel --prefix=$HOME --with-zlib --disable-shared --enable-fortran 
  #NOTE: This is with Intel MPI, if you have to use a different MPI, please adjust mpiicc and mpiifort accordingly to the scripts that call icc and ifort!
  #MPICH note: CC=mpicc FC=mpif90 ./configure --enable-parallel --prefix=$HOME --with-zlib --disable-shared --enable-fortran
  make
  make install
  cd ..

Installing netCDF
-----------------

.. code-block:: bash

  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz
  tar -xaf netcdf-4.3.0.tar.gz
  cd netcdf-4.3.0
  CPPFLAGS=-I$HOME/include LDFLAGS=-L$HOME/lib CC=mpiicc ./configure --enable-shared=no --prefix=$HOME 
  #MPICH note: CPPFLAGS=-I$HOME/include LDFLAGS=-L$HOME/lib CC=mpicc ./configure --enable-shared=no --prefix=$HOME 
  #NOTE: Check for this line to make sure netCDF is build with parallel I/O: 
  #"checking whether parallel I/O features are to be included... yes" This line comes at the very end (last 50 lines of configure run)!
  make
  make install
  cd ..


Installing Libxsmm
------------------

.. code-block:: bash

   git clone https://github.com/hfp/libxsmm
   cd libxsmm
   make generator
   #add libxsmm/bin to path


Installing Metis (Optional:PUML mesh format)
--------------------------------------------

.. code-block:: bash

   tar -xvf metis-5.1.0.tar.gz
   cd metis-5.1.0
   #edit ./include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
   make config cc=mpiicc cxx=mpiicpc prefix=$HOME #(must confirm compiler!!!)
   make install
   cd ..


Installing ParMetis (Optional:PUML mesh format)
-----------------------------------------------

.. code-block:: bash

  tar -xvf parmetis-4.0.3
  cd parmetis-4.0.3
  #edit ./metis/include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
  make config cc=mpiicc cxx=mpiicpc prefix=$HOME 
  #must confirm compiler! Better use the same path as metis!)
  make install
  cd ..

(Make sure parmetis/include contains metis.h and parmetis/lib contains
libmetis.a. Ortherwise, compile error: cannot find parmetis.)


Installing ASAGI (Optional)
---------------------------

See section :ref:`Installing ASAGI <installing_ASAGI>`.

.. _compiling-seissol:

Compiling SeisSol
-----------------

Get the latest version of SeisSol on git by cloning the whole repository
including all submodules:

.. code-block:: bash

   git clone https://github.com/SeisSol/SeisSol.git
   git submodule update --init

Add the following build variables to the file
build/options/supermuc_mac_cluster.py

.. code-block:: python

   compileMode='release' 
   parallelization='hybrid' 
   arch='$ARCH' 
   order='$ORDER' 
   generatedKernels = 'yes'
   compiler = 'intel'
   logLevel = 'info'

   netcdf='yes' 
   netcdfDir='path_to_netcdf' 
   hdf5='yes'
   hdf5Dir='path_to_hdf5'

   ##  additionally for puml mesh format
   metis = 'yes'
   metisDir='path_to_parmetis'

   ##  optional for ASAGI
   asagi = 'yes'
   zlibDir = 'path_to_asagi' #e.g. <path_to_ASAGI>/build/lib/

| with: 
| compileMode - release / relWithDebInfo/ debug
| parallelization - omp/ mpi / hybrid (mpi/openmp)
| logLevel - info/ debug, warning or error 
| ARCH - target architecture 
| ORDER - convergence order (=max polynomial order +1)
| generatedKernels - yes/no

Get your executable with

.. code-block:: bash

   scons -j 32 buildVariablesFile=build/options/supermuc_mac_cluster.py

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
   ``./SeisSol_release_.... PARAMETER.PAR``

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also :ref:`A first example <a_first_example>`.
