.. _build_seissol:

Compiling and Running SeisSol
-----------------------------

Compiling SeisSol
~~~~~~~~~~~~~~~~~
Get the latest version of SeisSol on git by cloning the whole repository,
including all submodules:

.. code-block:: bash

   git clone --recursive https://github.com/SeisSol/SeisSol.git

Then, given that you have installed all required dependencies, you may compile SeisSol for CPUs as follows (or use ``ccmake`` for a GUI):

.. code-block:: bash

   mkdir build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=hsw -DPRECISION=double -DORDER=4 -DGEMM_TOOLS_LIST=Eigen -DEQUATIONS=elastic -DGRAPH_PARTITIONING_LIBS=none ..
   make -j 4

The following parameters are important to note:

``PRECISION`` is either ``single`` or ``double`` and corresponds to 32-bit floating point numbers and 64-bit floating point numbers
used during all calculations. If your simulation gives an error when using single precision, it can sometimes be fixed by re-compiling for double precision.
The ``ORDER`` parameter may range from 2 to 8. The ``EQUATIONS`` parameter may be ``elastic``, ``viscoelastic2``, ``poroelastic``,
or ``anisotropic``. For the ``viscoelastic2`` setting, you will need to add e.g. ``-DNUMBER_OF_MECHANISMS=3`` (or replace 3 by any other number greater than zero).

Once you have installed a GEMM generator, choose a different option in ``GEMM_TOOLS_LIST`` to obtain a faster version of SeisSol.
For example, ``GEMM_TOOLS_LIST=libxsmm`` makes libxsmm generate the kernels.
Also, you may set ``GRAPH_PARTITIONING_LIBS`` to a graph partitioning library: we currently allow the values
``parmetis``, ``parhip`` and ``ptscotch``.

The ``HOST_ARCH`` parameter tells the compiler for which CPUs to optimize for. Note in particular that using the wrong host architecture may lead to "illegal instruction" (SIGILL) errors.
The host architecture ``hsw`` (an abbreviation for Haswell) will generate at most AVX2 instructions and should usually work when installing SeisSol on a local Linux machine for testing purposes. On clusters or workstations using Intel Skylake CPUs or newer (e.g. SuperMUC), use `skx` instead.
For an overview, the following values are accepted by ``HOST_ARCH``:

.. list-table:: Title
   :widths: 20 40 40
   :header-rows: 1

   * - ``HOST_ARCH``
     - Architecture
     - Notes
   * - ``noarch``
     - No architecture-specific optimizations
     - Generates plain x86-64 instructions, without SIMD instructions like SSE/AVX/AMX etc.
   * - ``wsm``
     - Intel Westmere architecture
     - Generates SSE instructions (up to SSE 3).
   * - ``snb``
     - Intel Sandy Bridge architecture
     - Generates AVX instructions.
   * - ``hsw``
     - Intel Haswell
     - Generates AVX2 instructions.
   * - ``skx``
     - Intel Skylake-X (including Skylake-SP)
     - Generates AVX-512{F,CD,BW,DQ,VL} instructions. (NOTE: Skylake desktop processors are NOT included here, unless they contain an "X" in their name, such as e.g. i9 7800X)
   * - ``knc``
     - Intel Knight's Corner (Xeon Phi coprocessor)
     - Generates Knight's Corner-specific instructions.
   * - ``knl``
     - Intel Knight's Landing (Xeon Phi, optionally as coprocessor)
     - Generates AVX-512{F,CD,PF,ER} instructions.
   * - ``rome``
     - AMD Rome (Zen, 2nd generation)
     - Generates AVX2 instructions. For the libxsmm kernel generator, it is deemed equivalent to `hsw`.
   * - ``thunderx2t99``
     - ARM ThunderX2
     - 
   * - ``power9``
     - IBM PowerPC 9
     - 

In case of a manual installation of dependencies (and not using the `setup.sh` script from the "Installing Dependencies" page), you may have to prepend :code:`CMAKE_PREFIX_PATH` and :code:`PKG_CONFIG_PATH` to the cmake command, e.g. for dependencies installed in :code:`${PREFIX}`:

.. code-block:: bash

    CMAKE_PREFIX_PATH=$PREFIX:$CMAKE_PREFIX_PATH PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig/:$PKG_CONFIG_PATH CC=...

It is also important that the executables of the matrix multiplication generators (Libxsmm, PSpaMM) have to be in your :code:`$PATH`.

You can also compile just the proxy by ``make SeisSol-proxy`` or only SeisSol with ``make SeisSol-bin``.

You can also run ``ccmake ..`` to see all available options and toggle them.

.. figure:: LatexFigures/ccmake.png
   :alt: An example of ccmake with some options



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
It is named :code:`./SeisSol_Release_....`. As argument, give it a SeisSol parameter file.

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also :ref:`A first example <a_first_example>`.
