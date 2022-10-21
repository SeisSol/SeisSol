Compiling and running SeisSol
-----------------------------

Compiling SeisSol
~~~~~~~~~~~~~~~~~
Get the latest version of SeisSol on git by cloning the whole repository
including all submodules:

.. code-block:: bash

   git clone --recursive https://github.com/SeisSol/SeisSol.git

If you have compiled :ref:`seissol-env with spack<spack_installation>`, load the module and compile SeisSol with (e.g.):

.. code-block:: bash

   mkdir build-release && cd build-release
   CC=mpicc CXX=mpiCC FC=mpif90 cmake -DCOMMTHREAD=ON -DNUMA_AWARE_PINNING=ON -DASAGI=ON -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=skx -DPRECISION=double -DORDER=4 -DGEMM_TOOLS_LIST=LIBXSMM,PSpaMM ..
   make -j 4

Please adapt ``CC``, ``CXX`` and ``FC`` to the mpi compilers you used for compiling the dependencies.
In case of a manual installation of dependencies, you may have to prepend :code:`CMAKE_PREFIX_PATH` and :code:`PKG_CONFIG_PATH` to the cmake command, e.g. for dependencies installed in :code:`${HOME}`:

.. code-block:: bash

    CMAKE_PREFIX_PATH=~:$CMAKE_PREFIX_PATH PKG_CONFIG_PATH=~/lib/pkgconfig/:$PKG_CONFIG_PATH CC=...

It is also important that the executables of the matrix multiplication generators (Libxsmm, PSpaMM) have to be in :code:`$PATH`.


You can also compile just the proxy by ``make SeisSol-proxy`` or only SeisSol with ``make SeisSol-bin`` 

Note: CMake tries to detect the correct MPI wrappers.

You can also run ``ccmake ..`` to see all available options and toggle them.

.. figure:: LatexFigures/ccmake.png
   :alt: An example of ccmake with some options

Compile with Score-P
""""""""""""""""""""

The Score-P measurement infrastructure is a highly scalable and easy-to-use tool suite for profiling and event tracing of HPC applications.
To compile with Score-P, use:

.. code-block:: bash

    SCOREP_WRAPPER=off CXX=scorep-mpic++ CC=scorep-mpicc FC=scorep-mpif90 cmake ..
    SCOREP_WRAPPER_INSTRUMENTER_FLAGS="--user --thread=omp --nomemory" make

Running SeisSol
~~~~~~~~~~~~~~~

1. Follow the instructions on :ref:`Configuration <Configuration>`.
2. Run SeisSol version of interest. To run the example:
   :code:`./SeisSol_Release_.... parameter.par`

Further information regarding meshing and parameter files etc. can be
found in the documentation folder. See also :ref:`A first example <a_first_example>`.
