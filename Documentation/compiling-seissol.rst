.. _build_seissol:

Compiling and Running SeisSol
-----------------------------

For this page, we will assume that all the environment variables
from the previous section to be loaded (cf. section ). That is, you have your SeisSol build directory in ``$SEISSOL_BASE``.

Get the latest version of SeisSol on git by cloning the whole repository,
including all submodules:

.. code-block:: bash

   git clone --recursive https://github.com/SeisSol/SeisSol.git

You may also do a shallow clone by adding ``--depth=1`` which will save a bit of download volume.

Building SeisSol
~~~~~~~~~~~~~~~~
Then, given that you have installed all required dependencies, you may compile SeisSol for CPUs as follows (or use ``ccmake`` for a GUI):

.. code-block:: bash

   mkdir build && cd build
   cmake -DNUMA_AWARE_PINNING=ON -DASAGI=ON -DPRECISION=double -DORDER=4 -DEQUATIONS=elastic ..
   make -j 4

Note that you may need to adjust your host architecture, especially if you are on an ARM-based machine: then you will add at least ``-DHOST_ARCH=neon``, or one of the other ARM-compatible architectures.

A list of build flags can be found here.

In case of a manual installation of dependencies (and not using the ``setup.sh`` script from the "Installing Dependencies" page), you may have to prepend :code:``CMAKE_PREFIX_PATH`` and :code:``PKG_CONFIG_PATH`` to the cmake command, e.g. for dependencies installed in :code:``${PREFIX}``:

.. code-block:: bash

    CMAKE_PREFIX_PATH=$SEISSOL_PREFIX:$CMAKE_PREFIX_PATH PKG_CONFIG_PATH=$SEISSOL_PREFIX/lib/pkgconfig/:$PKG_CONFIG_PATH

It is also important that the executables of the matrix multiplication generators (Libxsmm, PSpaMM) have to be in your :code:``$PATH``.

You can also compile just the proxy by ``make SeisSol-proxy`` or only SeisSol with ``make SeisSol-bin``.

You can also run ``ccmake ..`` to see all available options and toggle them.

.. figure:: LatexFigures/ccmake.png
   :alt: An example of ccmake with some options

Building for GPUs
~~~~~~~~~~~~~~~~~

When building for GPUs, you will need to install SYCL and gemmforge, chainforge as well.
Also, you will need to supply the device backend and the device architecture.
The backend can be ``cuda``, ``hip``, or one of the SYCL implementations, ``hipsycl`` or ``oneapi``.

Generally, we recommend compiling for simple precision on GPUs.

To give some examples: for an NVIDIA GH200 Superchip, you would therefore need to set (note that the ``HOST_ARCH`` field here needs to be changed to something ARM-based)

.. code-block:: bash

   mkdir build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=sve128 -DPRECISION=single -DORDER=4 -DDEVICE_BACKEND=cuda -DDEVICE_ARCH=sm_90 ..
   make -j 4

For an AMD Instinct MI250X GPU with Zen 3 CPU, you could do

.. code-block:: bash

   mkdir build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=milan -DPRECISION=single -DORDER=4 -DDEVICE_BACKEND=hip -DDEVICE_ARCH=gfx90a ..
   make -j 4

On an Intel Data Center Max GPU 1550, you could set

.. code-block:: bash

   mkdir build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=skx -DPRECISION=single -DORDER=4 -DDEVICE_BACKEND=oneapi -DSYCL_CC=dpcpp -DDEVICE_ARCH=pvc ..
   make -j 4

Cray Compiler Environments and GPUs
"""""""""""""""""""""""""""""""""""

When compiling AdaptiveCpp with the Cray Compiler Environment, it may not find some MPI files.
Therefore, you can use the following workaround (given that you know the base compilers, here written as ``COMPILER_C`` and  ``COMPILER_CXX``):

.. code-block:: bash

   mkdir build && cd build
   CC=$COMPILER_C CXX=$COMPILER_CXX CFLAGS=$(cc --cray-print-opts=all) CXXFLAGS=$(CC --cray-print-opts=all) cmake $!
   make -j 4

Why different builds?
~~~~~~~~~~~~~~~~~~~~~

Currently, SeisSol builds have the following constraints: they are restricted to one PDE, one precision and one polynomial degree usage for discretization

* a single equation system (isotropic elastic, anisotropic elastic, viscoelastic, poroelastic)
* a single polynomial discretization degree (2 to 7)
* a precision (float or double)
* a target architecture

Subsequently, it can be useful to re-build SeisSol multiple times with different configurations.
Each of these SeisSol builds has a different executable name, and they can be installed side-by-side.

Finding out your target architecture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For SeisSol to work optimally, you will need to find out your CPU and GPU architecture
you want to run on. That is, if you have a cluster, you will usually find the specifications
within the documentation of it. A list of the supported architectures can be found under TODO.

Generally speaking, if you encounter ``SIGILL`` errors, change your ``HOST_ARCH`` to a less demanding one (e.g. ``skx`` to ``hsw``).

A few heuristics may help in the beginning:

* if you work with your personal computer or laptop, you will be good with ``hsw`` or ``rome``.\footnote{If your computer is very old (i.e. 2013 or earlier), then you may have to check out snb, wsm or noarch instead} 
* on a cluster, or a high-end workstation (which supports AVX-512, or AVX10/512), you may use ``skx``.
* on an ARM machine, use ``neon`` and specify your CPU over the ``-mcpu`` parameter. If your machine supports SVE (such as A64FX or the Nvidia Grace CPU), then you can also use ``sve128``, ``sve256``, or ``sve512``; but you will still need to specify ``-mcpu``.
* if nothing else works, try out ``noarch``.

It shall be noted that support for the latest Apple Macbooks using Apple M1 or M2 processors is highly experimental and may lead to an unstable build or incorrect results.

For a list of known CPU configurations, check out TODO.


For GPUs, you may determine the local GPU if you have a viable ``llvm``/``clang`` installation available, e.g. by loading a module.
Then you can run the following commands.
If you are on a cluster, make sure to run them on a compute node.
* run ``nvptx-arch`` to determine the architecture version of the Nvidia GPUs
* run ``amdgpu-arch`` to determine the architecture version of the AMD GPUs (including the integrated GPUs of AMD CPUs)

Alternatively, you can also use:
* ``nvidia-smi --query-gpu compute_cap --format=csv`` for Nvidia GPUs. The numbers will be printed in the format "x.y" which corresponds to "sm_xy". E.g., "8.6" will become "sm_86".
* ``clinfo -l`` for AMD GPUs or Intel GPUs.
* ``rocminfo | grep gfx`` for AMD GPUs.

Once again, see TODO for a list of host arch parameters.

Compile with Score-P
""""""""""""""""""""

The Score-P measurement infrastructure is a highly scalable and easy-to-use tool suite for profiling and event tracing of HPC applications.
To compile with Score-P, use:

.. code-block:: bash

    SCOREP_WRAPPER=off CXX=scorep-mpic++ CC=scorep-mpicc cmake ..
    SCOREP_WRAPPER_INSTRUMENTER_FLAGS="--user --thread=omp --nomemory" make

Running SeisSol
~~~~~~~~~~~~~~~

Once SeisSol has been compiled successfully, enter your build directory and run the SeisSol version of choice.
It is named :code:``./SeisSol_Release_....``. As argument, give it a SeisSol parameter file.

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also :ref:``A first example <a_first_example>``.
