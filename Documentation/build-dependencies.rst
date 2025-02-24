..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _build_dependencies:

Installing dependencies
=======================

For compiling SeisSol, you will need the following dependencies during build:

- A C++17-capable compiler. The following work:

  - GCC (>= 9.0; tested: 13.2)
  - ICX (tested: 2024.2)
  - Clang (tested: 18, 19)
  - NVHPC (tested: 24.09; currently still slow!)
  - Cray CE (however: no commthread support; needs ``SEISSOL_COMMTHREAD=0``)
  - ICC 2021.9 (except v1.3.0)
- CMake (>= 3.20)
- Python (>= 3.9)

  - numpy (>= 1.12.0)
  - setuptools (>= 0.61.0)

(the CI currently verifies the build against Gcc 13.2, Clang 19, ICX 2024.2, and NVHPC 24.09)

Additionally, you need the following libraries:

- MPI (Support for MPI Standard >= 2.2)
- Eigen (>= 3.4)
- Hdf5 (>= 1.8, parallel)
- easi (>= 1.5)
- (optional, recommended) a code generator

  - libxsmm (== 1.17 if using inline-assembly (LIBXSMM); otherwise >= 1.17 (LIBXSMM_JIT))
  - PSpaMM
- (optional, recommended) a mesh partitioning library

  - ParMETIS (needs METIS and sometimes GKLib to be installed as well)
  - ParHIP
  - PT-SCOTCH
- (optional, recommended) Netcdf (>= 4.4)
- (optional, recommended) ASAGI

All dependencies can be installed automatically with spack or manually.
For the maximal performance and functionality,
we always recommend installing all optional dependencies.
(they can be skipped when pursuing a purely minimal installation)

For the GPU version, the following packages need to be installed as well:

- A GPU runtime model.

  - CUDA (>= 11.8) for NVIDIA GPUs
  - ROCm with HIP (>= 6.0.0) for AMD GPUs
  - SYCL: either AdaptiveCpp (>= 23.10) or DPC++; for Intel, AMD, NVIDIA, or other GPUs (note: SYCL was mandatory up to including SeisSol 1.3.0)
- gemmforge (>= 0.0.218, for NVIDIA, AMD and Intel GPUs)
- chainforge (>= 0.0.3, for NVIDIA and AMD GPUs)

.. _spack_installation:

Installation with Spack
-----------------------

The `Spack <https://github.com/spack/spack/wiki>`_ repository contains ``seissol`` as a package which automatically installs all necessary dependencies.
All relevant :ref:`build parameters <build_parameters>` are mapped to Spack variants; and architecture information is mostly inferred automatically—as known by Spack.

In case you work on a cluster with an older Spack version which does not have SeisSol as a package yet,
we also provide an out-of-tree dependency setup environment under https://github.com/SeisSol/seissol-spack-aid/tree/main/spack as ``seissol-env``.
See for reference our documentation on how to compile seissol-env on :ref:`SuperMUC-NG <compile_run_supermuc>`, :ref:`Shaheen <compile_run_shaheen>` (Cray system) and :ref:`Frontera <compile_run_frontera>`.
However, ``seissol-env`` is deprecated; we strongly advise switching to the upstream Spack package instead.

Manual installation
-------------------

We recommend setting up a folder on your system which will contain the build files for all dependencies.
For example, do ``mkdir ~/seissol; cd seissol``.

In _all_ cases before building something manually,
make sure to check the installed module files. That is, type ``module avail`` and look for the software you want to use.
Type ``module load NAME`` to load the respective software (with ``NAME`` being the name of the software, including the text after the slash, e.g. ``module load cmake/3.20.0``).

However, in some cases, the modules may be incomplete. Check that especially when using NVHPC, or the components for building AdaptiveCpp (LLVM, Boost).

.. _build_env:

Setting helpful environment variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a file ``setup.sh`` with the following enviroment variables. The following script assumes that you use the folder `~/seissol`
for building—adjust to your actual location.

.. code-block:: bash

  # For the Intel compiler
  # source /opt/intel/compiler/VERSION/bin/compilervars.sh intel64

  # write the path here which you created your directory in (you can figure it out via the `pwd` command)
  # here, $HOME/my-seissol-installation is used as an example; customize to your likening
  export SEISSOL_BASE=$HOME/my-seissol-installation

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

  # run "source ~/my-seissol-installation/setup.sh" to apply environment to the current shell

Required dependencies
~~~~~~~~~~~~~~~~~~~~~

We assume that you have a compiler already installed. The same goes for a suitable Python installation.
You will also need CMake in version 3.20.0 or above. Most likely, you system will already have a
version of CMake installed; however, you may have to load a module to get a new enough version.

If you do not have CMake in a new enough version available, you may also install it manually as follows.

.. code-block:: bash

  (cd $(mktemp -d) && wget -qO- https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0-Linux-x86_64.tar.gz | tar -xvz -C "." && mv "./cmake-3.20.0-linux-x86_64" "${SEISSOL_PREFIX}/bin/cmake")

Note that this extracts CMake to the directory ``${SEISSOL_PREFIX}/bin/cmake``, if you wish you can adjust that path. Note that you may now also use ``ccmake`` to get a terminal UI for configuring the following libraries.

Required libraries
~~~~~~~~~~~~~~~~~~

The following libraries need to be installed for all SeisSol CPU and GPU builds.
To get a working CPU build, installing all libraries described here is enough.
However, installing a GEMM generator and a graph partitioner is still recommended for better performance and better load balancing, respectively.

Installing HDF5
"""""""""""""""

If your system does not have it e.g. as a module file (type ``module avail | grep hdf5`` to look for it),
you may compile it manually with the following commands:

.. code-block:: bash

  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.3/src/hdf5-1.12.3.tar.bz2
  tar -xaf hdf5-1.12.3.tar.bz2
  cd hdf5-1.12.3
  CPPFLAGS="-fPIC ${CPPFLAGS}" CC=mpicc CXX=mpicxx ./configure --enable-parallel --prefix=$SEISSOL_PREFIX --with-zlib --disable-shared
  make -j8
  make install
  cd ..

Make sure to use the MPI compiler wrappers here. For the Intel compilers, you may use ``CC=mpiicx CXX=mpiicpx`` instead.

HDF5 is used for both mesh input (the PUML format, default in SeisSol) and high-order mesh output, as well as for checkpointing.

Installing Eigen
""""""""""""""""

Uf you do not have Eigen installed, you may do so manually as follows:

.. code-block:: bash

   wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
   tar -xf eigen-3.4.0.tar.gz
   cd eigen-3.4.0
   mkdir build && cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
   make install
   cd ../..

Eigen conveniently uses CMake as a build system for itself.
It is used in SeisSol for setting up matrices and other numerical computations, and optionally, also as code generator backend for matrix chain products.

Installing Easi
"""""""""""""""

Easi is used for setting up the model parameters.
It is (most likely) not already installed on your system or as a module file, as it is a more SeisSol-specific library.
You can find the installation instructions for it `under this link <https://easyinit.readthedocs.io/en/latest/getting_started.html>`_.

And with that, we're good to go!

Code generators for CPUs (optional, recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We support the following CPU code generators:

- libxsmm (``libxsmm\_gemm\_generator``) will give reasonable performance on most ``x86`` machines. Its JIT variant also supports ARM CPUs.
- PSpaMM (``pspamm-generator``): can handle some special cases faster; recommended mostly on AVX512-capable machines in conjunction with LIBXSMM. Otherwise slightly slower than LIBXSMM.
- Eigen: should work on all available architectures, but slower. Recommended, if you have trouble with the afore-mentioned code generators.

Note that using Eigen does not result in any additional dependencies, since it is needed in SeisSol anyways.

These GEMM generators are used to create optimized code for small matrix-matrix multiplications; as such their requirements differ from the usually-used BLAS libraries.

For GPU code generators, we currently only support gemmforge and chainforge, and the latter (chainforge) is recommended.

Installing Libxsmm (CPU)
""""""""""""""""""""""""

.. code-block:: bash

   git clone --depth=1 --branch 1.17 https://github.com/libxsmm/libxsmm
   cd libxsmm
   make generator
   cp bin/libxsmm_gemm_generator $SEISSOL_PREFIX/bin/
   cd ..

Note that you need to use version 1.17; newer versions will not work with SeisSol.

.. _installing_pspamm:

Installing PSpaMM (CPU)
"""""""""""""""""""""""

PSpaMM is a Python package, meaning that you can directly install it via pip:

.. code-block:: bash

   pip3 install --user git+https://github.com/SeisSol/PSpaMM.git

Usually PSpaMM is fast, but a bit slower than LIBXSMM. However, in some cases, it supersedes it.

Mesh partitioning library (optional, recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a good load balance on large clusters, SeisSol utilizes a mesh partitioning library during the startup of the simulation.
Currently, the software supports the following libraries:

-  ParMETIS (compile with ``IDXTYPEWIDTH=64``)
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
  sed -i 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g'  ./metis/include/metis.h
  make config cc=mpicc cxx=mpicxx prefix=$SEISSOL_PREFIX
  make install
  cp build/Linux-x86_64/libmetis/libmetis.a $SEISSOL_PREFIX/lib
  cp metis/include/metis.h $SEISSOL_PREFIX/include
  cd ..

Again, make sure to use the MPI compiler wrappers here, and adjust accordingly, as with the Hdf5 installation.

Also, make sure ``$SEISSOL_PREFIX/include`` contains ``metis.h`` and ``$SEISSOL_PREFIX/lib`` contains
``libmetis.a``. Otherwise, a compile error may come up.

Other functionalities (optional, recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

netCDF
""""""

NetCDF is needed for convergence tests, as these use periodic boundary conditions, and such are not yet supported by the PUML mesh format.
Also, point sources utilize the netCDF backend for one type of them.
Once again, if you do not have it installed (sometimes it comes bundled with HDF5), you may do so manually.

.. code-block:: bash

  wget https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/netcdf-c-4.8.1.tar.gz
  tar -xaf netcdf-4.8.1.tar.gz
  cd netcdf-4.8.1
  CFLAGS="-fPIC ${CFLAGS}" CC=h5pcc ./configure --enable-shared=no --prefix=$SEISSOL_PREFIX --disable-dap
  #NOTE: Check for this line to make sure netCDF is built with parallel I/O:
  #"checking whether parallel I/O features are to be included... yes" This line comes at the very end (last 50 lines of configure run)!
  make -j8
  make install
  cd ..

Note that for ``h5pcc`` to exist, you need to have compiled Hdf5 with MPI support. Using ``h5cc`` will most likely not work.

ASAGI
"""""

See section :ref:`Installing ASAGI <installing_ASAGI>`. A working parallel Netcdf installation is required for ASAGI.
Furthermore, you will need to compile easi with ASAGI support, by setting ``ASAGI=ON`` in the easi CMake after having installed ASAGI.

Additional requirements for GPUs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For GPUs, we need some more packages.

Installing GemmForge, ChainForge
""""""""""""""""""""""""""""""""

.. _gemmforge_installation:

The GPU code generators are called GemmForge and ChainForge.
Conveniently, they come as Python packages and can be installed with the following commands.

.. code-block:: bash

   pip3 install --user git+https://github.com/SeisSol/gemmforge.git
   pip3 install --user git+https://github.com/SeisSol/chainforge.git

Note that ChainForge is optional, but highly recommended for AMD and NVIDIA GPUs.
However, it does currently not support code generation to SYCL.

Once you have GemmForge/ChainForge ready, you are set for compiling SeisSol with GPUs.

Installing SYCL (for GPUs; optional for AMD and NVIDIA GPUs)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

See section :ref:`Installing SYCL <installing_SYCL>`.

SYCL is necessary for non-NVIDIA and non-AMD GPUs.
But you may (optionally) also compile SeisSol to use SYCL for NVIDIA or AMD GPUs.

Compiling SeisSol
-----------------

And with that, we're ready to compile SeisSol itself. For that, proceed to the next page
:ref:`Compiling SeisSol <build_seissol>`.
