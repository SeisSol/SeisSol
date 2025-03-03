..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _build_parameters:

Build parameter reference
=========================

**For Developers:** if you add a CMake parameter, make sure to document it here.

For a good overview and easy, we recommend the use of ``ccmake``. To do so, simply replace the use of ``cmake`` by ``ccmake``.
(i.e. add a "c" in front of the command)
Alternatively, you may also modify the ``CMakeCache.txt`` in your build directory directly (for it to exist, ``cmake`` needs to be run first).

Simulation- and optimization-specific parameters
------------------------------------------------

The following parameters will alter the name of the SeisSol executable.
You may explicitly compile and install multiple of these configurations at the same time—however, you will have to re-do the CMake and build process for each configuration.

- ``CMAKE_BUILD_TYPE``: Either ``Release``, ``RelWithDebInfo``, or ``Debug``

    * ``Release``: default value; usually optimizes with ``-O3``.
    * ``RelWithDebInfo``: optimizes a bit less than ``Release`` (usually ``-O2``), but offers more debug information.
    * ``Debug``: Also enables assertions. Note that more bugs can appear here than when running with the release options.
- ``EQUATION``: the equation system to compile for

    * ``elastic``: isotropic elastic
    * ``viscoelastic``: obsolete, use ``viscoelastic2`` instead
    * ``viscoelastic2``: isotropic viscoelastic. Requires a positive number of relaxation mechanisms (``NUMBER_OF_MECHANISMS``).
    * ``anisotropic``: anisotropic elastic; essentially uses the same kernels as the ``elastic`` mode, but with more general matrices
    * ``poroelastic``: isotropic poroelastic
- ``NUMBER_OF_MECHANISMS``: the number of mechanisms for viscoelastic simulations. For all other equations, this parameter is required to be 0.
- ``ORDER``: the expected convergence order. It corresponds to the polynomial degree plus 1. The order is used for both space and time integration. For example, if you specify order 4, you will be using polynomials of degree 3 in space and time. Note that a higher order can impact the performance greatly.
- ``PRECISION``:

    * ``single``: use single precision. Recommended in general for faster simulations. But especially for consumer GPUs (i.e. Nvidia GeForce, AMD Radeon, Intel ARC etc.), since these have usually a high performance difference between single and double precision.
    * ``double``: use double precision. Recommended, if your simulation fails with Inf/NaN errors in single precision builds. See also https://github.com/SeisSol/SeisSol/issues/200 .

Besides these, the host or, if enabled, the device architecture and backend are also encoded in the name of the executable.

Generic parameters
------------------

- ``LOG_LEVEL_MASTER``: determines the minimum relevance level of log messages that are printed. Possible values are ``error``, ``warning``, ``info``, and ``debug``.
- ``ASAGI``: enable or disable ASAGI as input for fault and material information. Note that easi will need to be built with ASAGI support if you select this option.
- ``COVERAGE``: determine the code coverage (currently only relevant for the CI).
- ``DR_QUAD_RULE``: the quadrature rule used for the Dynamic Rupture. Currently permits either ``dunavant`` or ``stroud``, referring to the respective quadrature scheme names.
- ``PLASTICITY_METHOD``: changes the plasticity matrices to be used. Options are

    * ``nb``: nodal basis
    * ``ib``: interpolating basis
- ``PROXY_PYBINDING``: compile Python bindings for the SeisSol proxy
- ``TESTING``: compile unit tests for SeisSol
- ``TESTING_GENERATED``: compile unit tests
- ``LTO``: enable link-time optimization
- ``GRAPH_PARTITIONING_LIBS``: compile for the given graph partitioning libraries. Allowed options are:

    * ``parmetis``: Repository: https://github.com/KarypisLab/ParMETIS
    * ``PTScotch``: Repository:
    * ``ParHIP``: Repository: https://github.com/KaHIP/KaHIP
    * ``none``: Do not require a graph partitioning library. Recommended only for single-node/proxy builds.

CPU-specific parameters
-----------------------

- ``HOST_ARCH``: the parameter to tune the architecture for. See build-archs for an overview. Applies both compiler options and GEMM generation. If you get an illegal instruction error, check this variable first. Since there are many options, we have an own page for them. In short, use ``hsw`` when using a Linux/Windows PC or an older Mac. Use ``apple-m1`` or ``apple-m2`` on the latest Mac computers.
- ``MEMKIND``: enables HBM support.
- ``NUMA_AWARE_PINNING``: pin the free CPUs (those used for the communication and IO threads) according to the given NUMA domains.
- ``MEMORY_LAYOUT``: the sparsity patterns to apply. If not given, it will be chosen by the CPU architecture.
- ``GEMM_TOOLS_LIST``: the list for CPU GEMM generators that are used. Note that SeisSol benefits from code generation specifically for small matrices. Currently supports combinations of the following:

    * ``auto``: automatically selects the installed and usable GEMM generators out of ``LIBXSMM_JIT``, ``LIBXSMM`` and ``PSpaMM`` (in this order).
    * ``LIBXSMM_JIT``: libxsmm, in JIT mode
    * ``LIBXSMM``: libxsmm, in code generation mode
    * ``PSpaMM``: PSpaMM (may support more matrix types and offers higher-performance sparsity compared to libxsmm, slightly slower for dense matrices at the time being).
    * ``MKL``: Intel MKL
    * ``OpenBlas``: OpenBLAS

GPU-specific parameters
-----------------------

- ``DEVICE_BACKEND``: enables or disables the GPU backend. The following backends are available at the moment:

    * ``cuda``: Nvidia CUDA
    * ``hip``: AMD HIP, using ROCm (or CUDA).
    * ``hipsycl``: SYCL, more specifically AdaptiveCpp, formerly known as Open SYCL and hipSYCL. Provides support for Intel, AMD, and Nvidia GPUs. Repository: https://github.com/AdaptiveCpp/AdaptiveCpp
    * ``oneapi``: SYCL, more specifically Intel Data Parallel C++ (DPC++). Provides support for Intel, AMD, and Nvidia GPUs. The open source variant is located under https://github.com/intel/llvm
- ``DEVICE_ARCH``: the parameter to tune and compile the kernels for. See build-archs for an overview.
- ``SYCLCC``: chooses the SYCL compiler used for the dynamic rupture and point source parts. Can be either AdaptiveCpp (``hipsycl``) or DPC++ (``dpcpp``); the description is the same as for the ``DEVICE_BACKEND``.
- ``SYCL_USE_NVHPC``: if AdaptiveCpp is compiled with NVHPC support, and we use NVHPC
- ``USE_GRAPH_CAPTURING``: if a compute graph feature is available, then use it. This is currently the case for CUDA (since 11.0) and HIP (requires ROCm 6.1 or higher). Compute graph support for SYCL is still experimental, although DPC++/oneAPI implements an extension for it
- ``ENABLE_PROFILING_MARKERS``: Currently available for CUDA and HIP

Options currently known to be broken
------------------------------------

The following options are available, but need to be left in the state that they are in. Not doing so will most likely break the build process or the software.

- ``INTEGRATE_QUANTITIES``: assumed to be always disabled. Currently broken; it will probably be replaced in some version soon—when we refactor the IO component of SeisSol.
- ``NUMBER_OF_FUSED_SIMULATIONS``: needs to be 0 or 1. Currently still broken for any higher number; but a fix is planned, cf. https://github.com/SeisSol/SeisSol/pull/385
