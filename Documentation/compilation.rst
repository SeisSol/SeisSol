Compilation
===========

In order to run SeisSol, you need to first install:

-  Python (>= 3.5)
-  Numpy (>= 1.12.0)
-  hdf5 (>= 1.8, for instructions see below)
-  netcdf (C-Release) (>= 4.4, for instructions see below)
-  Intel compiler (>= 18.0, icc, icpc, ifort) or GCC (>= 7.0, gcc, g++, gfortran)
-  Some MPI implementation (e.g. OpenMPI)
-  ParMETIS for partitioning
-  libxsmm (libxsmm\_gemm\_generator) for small matrix multiplications
-  PSpaMM (pspamm.py) for small sparse matrix multiplications (required only on Knights Landing or Skylake)
-  CMake (>3.15), for compiling submodules ImpalaJIT and yaml-cpp, and for SeisSol itself

Initial Adjustments to .bashrc
------------------------------

Add the following lines to your .bashrc (vi ~/.bashrc).

.. code-block:: bash

  # For intel compiler
  # source /opt/intel/compiler/VERSION/bin/compilervars.sh intel64
  
  export PATH=$HOME/bin:$PATH
  export LIBRARY_PATH=$HOME/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
  export PKG_CONFIG_PATH=$HOME/lib/pkgconfig:$PKG_CONFIG_PATH
  export CMAKE_PREFIX_PATH=$HOME
  export EDITOR=vi
  export CPATH=$HOME/include:$CPATH 

  # run "exec bash" or "source ~/.bashrc" to apply environment to the current shell

Installing CMake
----------------

.. code-block:: bash

  wget -qO- https://github.com/Kitware/CMake/releases/download/v3.16.4/cmake-3.16.4-Linux-x86_64.tar.gz | tar -xvz -C "/" && mv "/cmake-3.16.4-Linux-x86_64" "${HOME}/bin/cmake"
  ln -s ${HOME}/bin/cmake/bin/cmake ${HOME}/bin


Note that this extracts CMake to the directory ${HOME}/bin/cmake, if you wish you can adjust that path.
  
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

  wget https://syncandshare.lrz.de/dl/fiJNAokgbe2vNU66Ru17DAjT/netcdf-4.6.1.tar.gz
  tar -xaf netcdf-4.6.1.tar.gz
  cd netcdf-4.6.1
  CFLAGS="-fPIC ${CFLAGS}" CC=h5pcc ./configure --enable-shared=no --prefix=$HOME --disable-dap
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

.. _installing_pspamm:

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
libmetis.a. Otherwise, compile error: cannot find parmetis.)


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

First, you need to compile SeisSol's submodules:

.. code-block:: bash

   # Impala
   cd submodules/ImpalaJIT
   mkdir build && cd build
   CC=mpiicc CXX=mpiicpc FC=mpiifort cmake -DCMAKE_INSTALL_PREFIX=$HOME -- ..
   make -j12 install

   # yaml-cpp
   cd ../../yaml-cpp
   mkdir build && cd build
   CC=mpiicc CXX=mpiicpc FC=mpiifort cmake -DCMAKE_INSTALL_PREFIX=$HOME -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_CPP_BUILD_TESTS=OFF -- ..
   make -j12 install

   # go back to SeisSol root
   cd ../..


Compile SeisSol with (e.g.)

.. code-block:: bash

    mkdir build-release && cd build-release
    CC=mpiicc CXX=mpiicpc FC=mpiifort  CMAKE_PREFIX_PATH=~:$CMAKE_PREFIX_PATH PKG_CONFIG_PATH=~/lib/pkgconfig/:$PKG_CONFIG_PATH cmake -DNETCDF=ON -DMETIS=ON -DCOMMTHREAD=ON -DASAGI=OFF -DHDF5=ON -DCMAKE_BUILD_TYPE=Release -DTESTING=OFF  -DLOG_LEVEL=warning -DLOG_LEVEL_MASTER=info -DARCH=skx -DPRECISION=double ..
    make -j48

Here, the :code:`DCMAKE_INSTALL_PREFIX` controlls, in which folder the software is installed.
You have to adjust the :code:`CMAKE_PREFIX_PATH` and :code:`PKG_CONFIG_PATH` in the same manner - if you install all dependencies in a different directory, you need to replace :code:`${HOME}` by the path to this directory.
It is also important that the executables of the matrix mutiplication generators (Libxsmm, PSpaMM) have to be in :code:`$PATH`.
You can also compile just the proxy by :command:`make SeisSol-proxy` or only SeisSol with :command:`make SeisSol-bin`   

Note: CMake tries to detect the correct MPI wrappers.

You can also run :command:`ccmake ..` to see all available options and toggle them.

.. figure:: LatexFigures/ccmake.png
   :alt: An example of ccmake with some options


Running SeisSol
---------------

1. Follow the instructions on :ref:`Configuration <Configuration>`.
2. run SeisSol version of interest. To run the example:
   :command:`./SeisSol_release_.... PARAMETER.PAR`

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also :ref:`A first example <a_first_example>`.
