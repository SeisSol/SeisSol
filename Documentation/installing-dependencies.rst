Installing Dependencies
=======================

In order to run SeisSol, you need to first install:

-  Python (>= 3.5)
-  Numpy (>= 1.12.0)
-  hdf5 (>= 1.8, for instructions see below)
-  netcdf (C-Release) (>= 4.4, for instructions see below)
-  Intel compiler (>= 2021, icc, icpc, ifort) or GCC (>= 9.0, gcc, g++, gfortran)
-  Some MPI implementation (e.g. OpenMPI)
-  ParMETIS for partitioning (with IDXTYPEWIDTH=64)
-  libxsmm (libxsmm\_gemm\_generator) for small matrix multiplications
-  PSpaMM (pspamm.py) for small sparse matrix multiplications (required only on Knights Landing or Skylake)
-  CMake (>= 3.20) for the compilation of SeisSol

In addition, the following packages need to be installed for the GPU version of SeisSol:

- CUDA (>= 11.0)  for Nvidia GPUs, or HIP (ROCm>= 5.2.0) for AMD GPUs
- SYCL: either OpenSYCL >= 0.9.5 or DPC++
- gemmforge (>= 0.0.208)
- chainforge (>= 0.0.2, for Nvidia and AMD GPUs)


These dependencies can be installed automatically with spack or can be installed manually one by one.


.. _spack_installation:

Spack installation
------------------

`Spack <https://github.com/spack/spack/wiki>`_ is a HPC software package manager.
It automates the process of installing, upgrading, configuring, and removing computer programs.
In particular, our spack package `seissol-env` allows automatically installing all dependencies of SeisSol (e.g. mpi, hdf5, netcdf, easi, asagi, etc).
See https://github.com/SeisSol/seissol-spack-aid/tree/main/spack for details on the installation with spack.
See also for reference our documentation on how to compile seissol-env on :ref:`SuperMUC-NG <compile_run_supermuc>`, :ref:`Shaheen <compile_run_shaheen>` (Cray system) and :ref:`Frontera <compile_run_frontera>`.


Manual installation
-------------------

Initial Adjustments to .bashrc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~

.. code-block:: bash

  # you will need at least version 3.20.0 for GNU Compiler Collection 
  (cd $(mktemp -d) && wget -qO- https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0-Linux-x86_64.tar.gz | tar -xvz -C "." && mv "./cmake-3.20.0-linux-x86_64" "${HOME}/bin/cmake")
  
  # use version 3.16.2 for Intel Compiler Collection
  (cd $(mktemp -d) && wget -qO- https://github.com/Kitware/CMake/releases/download/v3.16.2/cmake-3.16.2-Linux-x86_64.tar.gz | tar -xvz -C "." && mv "./cmake-3.16.2-Linux-x86_64" "${HOME}/bin/cmake")
  
  ln -s ${HOME}/bin/cmake/bin/cmake ${HOME}/bin

Note that this extracts CMake to the directory ${HOME}/bin/cmake, if you wish you can adjust that path.
  
Installing HDF5
~~~~~~~~~~~~~~~

.. code-block:: bash

  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.8/src/hdf5-1.10.8.tar.bz2
  tar -xaf hdf5-1.10.8.tar.bz2
  cd hdf5-1.10.8
  CPPFLAGS="-fPIC ${CPPFLAGS}" CC=mpicc FC=mpif90 ./configure --enable-parallel --prefix=$HOME --with-zlib --disable-shared --enable-fortran 
  make -j8
  make install
  cd ..

Installing netCDF
~~~~~~~~~~~~~~~~~

.. code-block:: bash

  wget https://syncandshare.lrz.de/dl/fiJNAokgbe2vNU66Ru17DAjT/netcdf-4.6.1.tar.gz
  tar -xaf netcdf-4.6.1.tar.gz
  cd netcdf-4.6.1
  CFLAGS="-fPIC ${CFLAGS}" CC=h5pcc ./configure --enable-shared=no --prefix=$HOME --disable-dap
  #NOTE: Check for this line to make sure netCDF is built with parallel I/O: 
  #"checking whether parallel I/O features are to be included... yes" This line comes at the very end (last 50 lines of configure run)!
  make -j8
  make install
  cd ..

.. _installing_eigen3:

Installing Eigen3
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
   tar -xf eigen-3.4.0.tar.gz
   cd eigen-3.4.0
   mkdir build && cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=~
   make install
   cd ../..

.. _installing_libxsmm:

Installing Libxsmm
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone --branch 1.17 https://github.com/hfp/libxsmm
   cd libxsmm
   make generator
   cp bin/libxsmm_gemm_generator $HOME/bin
   cd ..

.. _installing_pspamm:

Installing PSpaMM
~~~~~~~~~~~~~~~~~



.. code-block:: bash

   git clone https://github.com/SeisSol/PSpaMM.git
   # make sure $HOME/bin exists or create it with "mkdir ~/bin"
   ln -s $(pwd)/PSpaMM/pspamm.py $HOME/bin/pspamm.py
   
Instead of linking, you could also add the following line to your .bashrc:

.. code-block:: bash

   export PATH=<Your_Path_to_PSpaMM>:$PATH


.. _installing_parmetis:

Installing ParMetis
~~~~~~~~~~~~~~~~~~~


.. code-block:: bash

  wget https://ftp.mcs.anl.gov/pub/pdetools/spack-pkgs/parmetis-4.0.3.tar.gz
  tar -xvf parmetis-4.0.3.tar.gz
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

See section :ref:`Installing ASAGI <installing_ASAGI>`.

.. _compiling-seissol:

Installing easi
~~~~~~~~~~~~~~~

Follow the `installation instructions <https://easyinit.readthedocs.io/en/latest/getting_started.html>`_.


Installing GemmForge, ChainForge (for GPUs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _gemmforge_installation:

.. code-block:: bash

   pip3 install --user git+https://github.com/ravil-mobile/gemmforge.git
   pip3 install --user git+https://github.com/ravil-mobile/chainforge.git

Installing SYCL (for GPUs)
~~~~~~~~~~~~~~~~~~~~~~~~~~

See section :ref:`Installing SYCL <installing_SYCL>`.


