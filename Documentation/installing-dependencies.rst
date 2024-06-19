Installing Dependencies
=======================

For compiling SeisSol, you will need the following dependencies during build:

- A C++17-capable compiler. For the 
  - GCC (>= 9.0)
  - ICC (>= 2021.0)
- CMake (>= 3.20)
- Python (>= 3.5)
- Numpy (>= 1.12.0)

Additionally, you need the following libraries:

- MPI (Support for MPI Standard >= 2.2)
- Eigen (>= 3.4)
- hdf5 (>= 1.8)
- easi (>= 1.2)
- (optional) a code generator

  - libxsmm (== 1.17, newer versions do not work with YATeTo right now)
  - PSpaMM
- (optional) mesh partitioning

  - ParMETIS (needs METIS to be installed as well)
  - ParHIP
  - PT-SCOTCH
- (optional) Netcdf (>= 4.4)
- (optional) ASAGI

All dependencies can be installed automatically with spack or manually by hand.
For a minimal installation,
you may avoid all optional dependencies. However, for the maximal performance, and for compute clusters,
we recommend having at least using a code generator and having a mesh partitioner linked to SeisSol.

For the GPU version, the following packages need to be installed as well:

- SYCL: either AdaptiveCpp (hipSYCL/Open SYCL) >= 0.9.3 or DPC++
- gemmforge (>= 0.0.207, for Nvidia, AMD and Intel GPUs)
- (optional) chainforge (>= 0.0.2, for Nvidia and AMD GPUs)
- (optional) CUDA (>= 11.0) for Nvidia GPUs, or HIP (ROCm>= 5.2.0) for AMD GPUs

.. _spack_installation:

Spack installation
------------------

The `Spack <https://github.com/spack/spack/wiki>`_ package `seissol-env` allows to automatically install all dependencies of SeisSol (e.g. mpi, hdf5, netcdf, easi, asagi, etc),
leaving only SeisSol itself left to be build.

See https://github.com/SeisSol/seissol-spack-aid/tree/main/spack for details on the installation with spack.
See also for reference our documentation on how to compile seissol-env on :ref:`SuperMUC-NG <compile_run_supermuc>`, :ref:`Shaheen <compile_run_shaheen>` (Cray system) and :ref:`Frontera <compile_run_frontera>`.

Manual Installation
-------------------

We recommend setting up a folder on your system which will contain the build files for all dependencies.
For example, do `mkdir ~/seissol; cd seissol`.

In _all_ cases before building something manually,
make sure to check the installed module files. That is, type `module avail` and look for the software you want to use.
Type `module load NAME` to load the respective software (with `NAME` being the name of the software, including the text after the slash, e.g. `module load cmake/3.20.0`).

Setting Helpful Environment Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a file `setup.sh` with the following enviroment variables. The following script assumes that you use the folder `~/seissol`
for building.

.. code-block:: bash

  # For intel compiler
  # source /opt/intel/compiler/VERSION/bin/compilervars.sh intel64
  
  # write the path here which you created your directory in (you can figure it out by typing `pwd`)
  export SEISSOL_BASE=$HOME/seissol

  export SEISSOL_PREFIX=$SEISSOL_BASE/local
  export PATH=$SEISSOL_PREFIX/bin:$PATH
  export LIBRARY_PATH=$SEISSOL_PREFIX/lib:$SEISSOL_PREFIX/lib64:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$SEISSOL_PREFIX/lib:$SEISSOL_PREFIX/lib64:$LD_LIBRARY_PATH
  export PKG_CONFIG_PATH=$SEISSOL_PREFIX/lib/pkgconfig:$SEISSOL_PREFIX/lib64/pkgconfig:$PKG_CONFIG_PATH
  export CMAKE_PREFIX_PATH=$SEISSOL_PREFIX:$CMAKE_PREFIX_PATH
  export CMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
  export CPATH=$SEISSOL_PREFIX/include:$CPATH 
  export C_INCLUDE_PATH=$SEISSOL_PREFIX/include:$C_INCLUDE_PATH
  export CXX_INCLUDE_PATH=$SEISSOL_PREFIX/include:$CXX_INCLUDE_PATH
  export EDITOR=nano # or e.g. vi,vim

  # run "source ~/seissol/setup.sh" to apply environment to the current shell

Required Dependencies
~~~~~~~~~~~~~~~~~~~~~

We assume that you have a compiler already installed. The same goes for a suitable Python installation.
You will also need CMake in version 3.20.0 or above. Most likely, you system will already have a
version of CMake installed; however, you may have to load a module to get a new enough version.

If you do not have CMake in a new enough version available, you may also install it manually as follows.

.. code-block:: bash

  (cd $(mktemp -d) && wget -qO- https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0-Linux-x86_64.tar.gz | tar -xvz -C "." && mv "./cmake-3.20.0-linux-x86_64" "${HOME}/bin/cmake")

Note that this extracts CMake to the directory ``${SEISSOL_PREFIX}/bin/cmake``, if you wish you can adjust that path. Note that you may now also use ``ccmake`` to get a terminal UI for configuring the following libraries.
  
Required Libraries
~~~~~~~~~~~~~~~~~~

The following libraries need to installed for all SeisSol CPU and GPU builds.
To get a working CPU build, installing all libraries described here is enough.
However, installing a GEMM generator and a graph partitioner is still recommended for better performance and better load balancing, respectively.

Installing HDF5
"""""""""""""""

We begin with HDF5 which is needed for reading meshes in PUML format (the default format in SeisSol, and the output of PUMgen).
If your system does not have it e.g. as a module file (type `module avail | grep hdf5` to look for it),
you may compile it manually with the following commands:

.. code-block:: bash

  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.8/src/hdf5-1.10.8.tar.bz2
  tar -xaf hdf5-1.10.8.tar.bz2
  cd hdf5-1.10.8
  CPPFLAGS="-fPIC ${CPPFLAGS}" CC=mpicc CXX=mpicxx ./configure --enable-parallel --prefix=$SEISSOL_PREFIX --with-zlib --disable-shared
  make -j8
  make install
  cd ..

Installing Eigen
""""""""""""""""

Next, we look at Eigen which conveniently uses CMake as a build system for itself.
Eigen is used in SeisSol for setting up matrices and other numerical computations, and optionally, also as code generator for matrix chain products.
Once again, if you do not have Eigen installed, you may do so manually as follows:

.. code-block:: bash

   wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
   tar -xf eigen-3.4.0.tar.gz
   cd eigen-3.4.0
   mkdir build && cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
   make install
   cd ../..

Lastly, we need easi which is (most likely) not already installed on your system or as a module file, as it is a more SeisSol-specific library.
It is used for setting up the model parameters.
Here you can find the `installation instructions <https://easyinit.readthedocs.io/en/latest/getting_started.html>`_.

And with that, we're good to go!

Additional Requirements for GPUs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For GPUs, we need some more packages.

Installing SYCL (for GPUs)
""""""""""""""""""""""""""

See section :ref:`Installing SYCL <installing_SYCL>`.

Installing GemmForge, ChainForge (for GPUs)
"""""""""""""""""""""""""""""""""""""""""""

.. _gemmforge_installation:

The GPU code generators are called GemmForge and ChainForge.
Conveniently, they come as Python packages and can be installed with the following commands.

.. code-block:: bash

   pip3 install --user git+https://github.com/SeisSol/gemmforge.git
   pip3 install --user git+https://github.com/SeisSol/chainforge.git

Note that ChainForge is optional.

Once you have SYCL and GemmForge (maybe also ChainForge) ready, you are set for compiling SeisSol with GPUs.

Code Generators for CPUs (optional, recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For CPU code generators, we support the following:

- Eigen
- libxsmm (libxsmm\_gemm\_generator) for small matrix multiplications
- PSpaMM (pspamm.py) for small sparse matrix multiplications (required only on Knights Landing or Skylake)

Note that using Eigen does not result in any additional dependencies, since it is needed in SeisSol anyways.
Other than that, we recommend using the combination libxsmm and PSpaMM.

For GPU code generators, we currently only support gemmforge and chainforge, and the latter (chainforge) is recommended.

Installing Libxsmm
""""""""""""""""""

(to save data, we only use a shallow clone)

.. code-block:: bash

   git clone --depth=1 --branch 1.17 https://github.com/hfp/libxsmm
   cd libxsmm
   make generator
   cp bin/libxsmm_gemm_generator $SEISSOL_PREFIX/bin/
   cd ..

.. _installing_pspamm:

Installing PSpaMM
"""""""""""""""""

(to save data, we only use a shallow clone)

.. code-block:: bash

   git clone --depth=1 https://github.com/SeisSol/PSpaMM.git
   # make sure $SEISSOL_PREFIX/bin exists or create it with "mkdir ~/bin"
   ln -s $(pwd)/PSpaMM/pspamm.py $SEISSOL_PREFIX/bin/pspamm.py

Mesh Partitioning (optional, recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a good load balance on large clusters, SeisSol utilizes a mesh partitioning library during the startup of the simulation.
Currently, the software supports the following libraries:

-  ParMETIS (compile with IDXTYPEWIDTH=64)
-  SCOTCH
-  ParHIP

The partitioning of SeisSol meshes with ParMETIS was tested in large simulations and is
generally recommended for academic users.
SCOTCH and ParHIP are free and open-source alternatives to ParMETIS and should be used by
users from industry or for-profit institutions (cf. `ParMETIS license <https://github.com/KarypisLab/ParMETIS/blob/main/LICENSE>`_).
A study comparing partition quality for SeisSol meshes can be found `here <https://home.in.tum.de/~schnelle/publications/bachelorsthesis-informatics-final.pdf>`_.

ParMETIS
""""""""

ParMETIS may be installed as follows:

.. code-block:: bash

  wget https://ftp.mcs.anl.gov/pub/pdetools/spack-pkgs/parmetis-4.0.3.tar.gz
  tar -xvf parmetis-4.0.3.tar.gz
  cd parmetis-4.0.3
  #edit ./metis/include/metis.h IDXTYPEWIDTH to be 64 (default is 32).
  make config cc=mpicc cxx=mpiCC prefix=$SEISSOL_PREFIX 
  make install
  cp build/Linux-x86_64/libmetis/libmetis.a $SEISSOL_PREFIX/lib
  cp metis/include/metis.h $SEISSOL_PREFIX/include
  cd ..

(Make sure $HOME/include contains metis.h and $HOME/lib contains
libmetis.a. Otherwise, compile error: cannot find parmetis.)

Other Software (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~

netCDF
""""""

NetCDF is needed for convergence tests, as these use periodic boundary conditions—and such are not supported by the PUML mesh format.
Also, point sources utilize the netCDF backend.
Once again, if you do not have it installed (sometimes it comes bundled with HDF5), you may do so manually.

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

ASAGI
"""""

See section :ref:`Installing ASAGI <installing_ASAGI>`.

Compiling SeisSol
-----------------

And with that, we're ready to compile SeisSol itself! For that, proceed to the next page
:ref:`Compiling SeisSol <build_seissol>`.
