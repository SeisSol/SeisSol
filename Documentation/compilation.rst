Compilation
===========

In order to run SeisSol, you need to first install:

-  Python (>= 3.5)
-  Numpy (>= 1.12.0)
-  SCons (>= 3.0, for instructions see below)
-  hdf5 (>= 1.8, for instructions see below)
-  netcdf (C-Release) (>= 4.4, for instructions see below)
-  Intel compiler (>= 17.0, icc, icpc, ifort) or GCC (>= 5.0, gcc, g++, gfortran)
-  Some MPI implementation (e.g. OpenMPI)
-  ParMETIS for partitioning
-  libxsmm (libxsmm\_gemm\_generator) for small matrix multiplications
-  PSpaMM (pspamm.py) for small sparse matrix multiplications (required only on Knights Landing or Skylake)
-  CMake (for compiling submodules ImpalaJIT and yaml-cpp)

Inital Adjustments to .bashrc
-----------------------------

Add the following lines to your .bashrc (vi ~/.bashrc).

.. code-block:: bash

  # For intel compiler
  # source /opt/intel/compiler/VERSION/bin/compilervars.sh intel64
  
  export PATH=$HOME/bin:$PATH
  export LIBRARY_PATH=$HOME/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
  export PKG_CONFIG_PATH=$HOME/lib/pkgconfig:$PKG_CONFIG_PATH
  export EDITOR=vi
  export CPATH=$HOME/include:$CPATH 

  # run "exec bash" or "source ~/.bashrc" to apply environment to current shell

Installing SCons
----------------

.. code-block:: bash

  wget http://prdownloads.sourceforge.net/scons/scons-3.0.5.tar.gz
  tar -xaf scons-3.0.5.tar.gz
  cd scons-3.0.5
  python setup.py install --prefix=$HOME
  cd ..


Installing HDF5
---------------

.. code-block:: bash

  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/hdf5-1.8.21.tar.bz2
  tar -xaf hdf5-1.8.21.tar.bz2
  cd hdf5-1.8.21
  CPPFLAGS="-fPIC ${CPPFLAGS}" CC=mpicc FC=mpif90 ./configure --enable-parallel --prefix=$HOME --with-zlib --disable-shared --enable-fortran 
  make -j8
  make install
  cd ..

Installing netCDF
-----------------

.. code-block:: bash

  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.6.1.tar.gz
  tar -xaf netcdf-4.6.1.tar.gz
  cd netcdf-4.6.1
  CFLAGS="-fPIC ${CFLAGS}" CC=h5pcc ./configure --enable-shared=no --prefix=$HOME 
  #NOTE: Check for this line to make sure netCDF is build with parallel I/O: 
  #"checking whether parallel I/O features are to be included... yes" This line comes at the very end (last 50 lines of configure run)!
  make -j8
  make install
  cd ..

.. _installing_libxsmm:

Installing Libxsmm
------------------

.. code-block:: bash

   git clone https://github.com/hfp/libxsmm
   cd libxsmm
   make generator
   cp bin/libxsmm_gemm_generator $HOME/bin
   cd ..

Installing PSpaMM
-----------------

.. code-block:: bash

   git clone https://github.com/peterwauligmann/PSpaMM.git
   ln -s $(pwd)/PSpaMM/pspamm.py $HOME/bin

Installing ParMetis (Optional: PUML mesh format)
------------------------------------------------

.. code-block:: bash

  wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
  tar -xvf parmetis-4.0.3
  cd parmetis-4.0.3
  #edit ./metis/include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
  make config cc=mpicc cxx=mpiCC prefix=$HOME 
  make install
  cp build/Linux-x86_64/libmetis/libmetis.a $HOME/lib
  cp metis/include/metis.h $HOME/include
  cd ..

(Make sure $HOME/include contains metis.h and $HOME/lib contains
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

Compile SeisSol with (e.g.)

.. code-block:: bash

  scons compiler=gcc netcdf=yes hdf5=yes order=4 parallelization=hybrid 

You may also save your favorite settings in a configuration file:
Add the following build variables to the file
build/options/supermuc_mac_cluster.py

.. code-block:: python

   compileMode='release' 
   parallelization='hybrid' 
   arch='$ARCH' 
   order='$ORDER' 
   generatedKernels = 'yes'
   compiler = 'gcc' # alternative: 'intel'
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

Running SeisSol
---------------

1. Follow the instructions on :ref:`Configuration <Configuration>`.
2. run SeisSol version of interest. To run the example:
   ``./SeisSol_release_.... PARAMETER.PAR``

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also :ref:`A first example <a_first_example>`.
