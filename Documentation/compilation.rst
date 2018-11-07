Compilation
===========

In order to run SeisSol, you need to first install: - Python (use
version 2.7 and above) - SCons (for instructions see below) - hdf5 (for
instructions see below) - netcdf (C-Release) (for instructions see
below) - Intel 2017 Compiler (icc/icpc/ifort) Suite (GCC support can be
added manually) - Some MPI implementation - ParMetis/Metis for
partitioning - Libxsmm package for small matrix multiplications - cmake
for compiling ImpalaJIT

Inital Adjustments to .bashrc
-----------------------------

1. add following lines to your .bashrc (vi ~/.bashrc), which should be
   on the NAS

   ::

       # These are the versions we validated against
       source /opt/intel/compiler/VERSION/bin/compilervars.sh intel64

       # User specific aliases and functions
       export PATH=$HOME/bin:$PATH
       export LD_LIBRARY_PATH=$HOME/lib:$HOME/lib64:$LD_LIBRARY_PATH
       export EDITOR=vi
       export SCONS_LIB_DIR=$HOME/lib64/scons-2.2.0/
       export PATH=$HOME/../libxsmm/bin:$PATH
       #This 2 lines have been suggested by @yzs981130 (not sure they are really necessary)
       export C_INCLUDE_PATH=$HOME/include:$C_INCLUDE_PATH 
       export CPLUS_INCLUDE_PATH=$HOME/include:$CPLUS_INCLUDE_PATH

       ######  ParMetis library necessary (Optional) ##############
       export PARMETIS_BASE='path_to_parmetis'
       export PARMETIS_LIBDIR='path_to_parmetis/lib'

Installing SCons
----------------

0. (wget http://prdownloads.sourceforge.net/scons/scons-2.2.0.tar.gz)
1. tar -xaf scons-2.2.0.tar.gz
2. cd scons-2.2.0
3. python setup.py install --prefix=$HOME
4. cd ..

Installing HDF5
---------------

0. (wget
   https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.11/src/hdf5-1.8.11.tar.bz2)
1. tar -xaf hdf5-1.8.11.tar.bz2
2. cd hdf5-1.8.11
3. CC=mpiicc FC=mpiifort ./configure --enable-parallel
   --prefix=:math:`HOME --with-zlib --disable-shared --enable-fortran    NOTE: This is with Intel MPI, if you have to use a different MPI, please adjust mpiicc and mpiifort accordingly to          the scripts that call icc and ifort!    MPICH note: CC=mpicc FC=mpif90 ./configure --enable-parallel --prefix=`\ HOME
   --with-zlib --disable-shared --enable-fortran
4. make
5. make install
6. cd ..

Installing netCDF
-----------------

0. (wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz)
1. tar -xaf netcdf-4.3.0.tar.gz
2. cd netcdf-4.3.0
3. CPPFLAGS=-I\ :math:`HOME/include LDFLAGS=-L`\ HOME/lib CC=mpiicc
   ./configure --enable-shared=no
   --prefix=:math:`HOME    MPICH note: CPPFLAGS=-I`\ HOME/include
   LDFLAGS=-L\ :math:`HOME/lib CC=mpicc ./configure --enable-shared=no --prefix=`\ HOME
   NOTE: Check for this line to make sure netCDF is build with parallel
   I/O: "checking whether parallel I/O features are to be included...
   yes" This line comes at the very end (last 50 lines of configure
   run)!
4. make
5. make install
6. cd ..

Installing Libxsmm
------------------

0. git clone https://github.com/hfp/libxsmm
1. cd libxsmm
2. make generator
3. add libxsmm/bin to path

Installing Metis (Optional:PUML mesh format)
--------------------------------------------

1. tar -xvf metis-5.1.0.tar.gz
2. cd metis-5.1.0
3. edit ./include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
4. make config cc=mpiicc cxx=mpiicpc prefix=$HOME (must confirm
   compiler!!!)
5. make install
6. cd ..

Installing ParMetis (Optional:PUML mesh format)
-----------------------------------------------

1. tar -xvf parmetis-4.0.3
2. cd parmetis-4.0.3
3. edit ./metis/include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
4. make config cc=mpiicc cxx=mpiicpc prefix=$HOME (must confirm
   compiler! Better use the same path as metis!)
5. make install
6. cd ..

(Make sure parmetis/include contains metis.h and parmetis/lib contains
libmetis.a. Ortherwise, compile error: cannot find parmetis.)

Compiling SeisSol
-----------------

Get the latest version of SeisSol on git by cloning the whole repository
including all submodules:

::

    git clone --recursive https://github.com/SeisSol/SeisSol.git

(Don't forget --recursive, otherwise the submodules are not cloned!)

Add the following build variables to the file
build/options/supermuc\_mac\_cluster.py

::

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
    hdf5Dir='path_to_netcdf'

    # additionally for puml mesh format
    metis = 'yes'
    metisDir='path_to_parmetis'

    # optional for ASAGI 
    zlibDir = 'path_to_ASAGI/lib'

with: compileMode - release / debug; parallelization - mpi / hybrid
(mpi/openmp); logLevel - info/ debug, warning or error; ARCH - target
architecture; ORDER - convergence order you want to use;
generatedKernels - yes/no; netcdfDir - path to netcdf; hdf5Dir - path to
hdf5; metisDir - path to parmetis; zlibDir - path to ASAGI lib
(optional) Get your executable with

::

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

1. Follow the instructions on:
   https://github.com/SeisSol/SeisSol/wiki/Configuration
2. run SeisSol version of interest 2a. To run the example:
   ./SeisSol\_release\_.... PARAMETER.PAR

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also [[A first example]].
