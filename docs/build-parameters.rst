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
- ``EQUATIONS``: the equation system to compile for

    * ``elastic``: isotropic elastic
    * ``acoustic``: acoustic
    * ``viscoelastic``: isotropic viscoelastic. Requires a positive number of relaxation mechanisms (``NUMBER_OF_MECHANISMS``). ``viscoelastic2`` is accepted as an alias.
    * ``viscoacoustic``: viscoacoustic. Requires a positive number of relaxation mechanisms (``NUMBER_OF_MECHANISMS``).
    * ``anisotropic``: anisotropic elastic; essentially uses the same kernels as the ``elastic`` mode, but with more general matrices
    * ``poroelastic``: isotropic poroelastic

- ``SOLVER`` (advanced): the scheme which advances the cells in time. It holds for every configuration of the executable whose material can be advanced with it; the others, and all with the default ``auto``, take the first scheme their equation system supports. So all configurations of one material in an executable share their scheme, and a configuration name does not carry it.

    * ``linearck``: the Cauchy-Kovalevskaya procedure for linear equations. The only choice for ``elastic``, ``acoustic`` and ``anisotropic``. For ``viscoelastic`` and ``viscoacoustic``, it carries the memory variables as additional quantities; this replaces ``VISCO_MODE=extend``.
    * ``linearckanelastic``: the Cauchy-Kovalevskaya procedure with the memory variables kept separately, which is faster. The default for ``viscoelastic`` and ``viscoacoustic``; this replaces ``VISCO_MODE=split``. On GPUs, it does not support free surfaces with gravity, Dirichlet or analytical boundaries; SeisSol stops during the initialization if the mesh has any.
    * ``stp``: the space-time predictor, which treats stiff source terms implicitly. The only choice for ``poroelastic``. It does not support free surfaces with gravity; SeisSol stops during the initialization if the mesh has any.

- ``NUMBER_OF_MECHANISMS``: the number of mechanisms for viscoelastic and viscoacoustic simulations. For all other equations, this parameter is required to be 0.
- ``ORDER``: the expected convergence order. It corresponds to the polynomial degree plus 1. The order is used for both space and time integration. For example, if you specify order 4, you will be using polynomials of degree 3 in space and time. Note that a higher order can impact the performance greatly.
- ``PRECISION``:

    * ``single``: use single precision. Recommended in general for faster simulations. But especially for consumer GPUs (i.e. Nvidia GeForce, AMD Radeon, Intel ARC etc.), since these have usually a high performance difference between single and double precision.
    * ``double``: use double precision. Recommended, if your simulation fails with Inf/NaN errors in single precision builds. See also https://github.com/SeisSol/SeisSol/issues/200 .

- ``NUMBER_OF_FUSED_SIMULATIONS``: the number of simulations run simultaneously.
- ``EXTRA_CONFIGS``: further configurations to build into the same executable, next to the one the parameters above describe, as a list of their names or of named sets of them (see ``CONFIGS``), e.g. ``elastic-o4-f32``. A name reads ``<material>[-m<mechanisms>]-o<order>-<f32|f64>[-s<fused simulations>]``, e.g. ``viscoelastic-m3-o6-f64`` or ``elastic-o3-f32-s8``; the mechanisms are left out when there are none, the fused simulations when there is only one. The scheme that advances a material (``SOLVER``) and the quadrature rule of the dynamic rupture faces (``DR_QUAD_RULE``) are the same for all configurations of an executable, so a name carries neither. A run chooses its configuration with ``Configuration`` in the ``equations`` section of the parameter file, or else with the environment variable ``SEISSOL_CONFIGURATION`` (see :ref:`env_vars`), and takes the first one otherwise; with ``ConfigMap``, mesh groups compute in configurations of their own (see :ref:`several_configurations`). All configurations share the alignment and the vector size of the host architecture.
- ``CONFIGS``: the configurations to build into the executable, as a list of their names, of named sets of them, or both, e.g. ``elastic-o4-f64;elastic-o4-f32`` or ``o4``. The named sets are in ``cmake/confignames.cmake``: ``o2`` to ``o8`` hold one order each, with the six materials in single and double precision, viscoelastic and viscoacoustic with 3 mechanisms; ``all`` holds these sets together (84 configurations), order 6 first; ``ci-cpu`` and ``ci-gpu`` hold the configurations the CPU and the GPU CI build, ``ci-sanitize`` those the sanitizer CI builds (the double precision half of ``ci-cpu``), and ``ci-interaction`` configurations of several orders and both precisions whose cells can be face neighbors (see :ref:`several_configurations`). The first configuration is the one a run takes if neither its parameter file nor ``SEISSOL_CONFIGURATION`` names one. If given, it replaces ``EQUATIONS``, ``NUMBER_OF_MECHANISMS``, ``ORDER``, ``PRECISION``, ``NUMBER_OF_FUSED_SIMULATIONS`` and ``EXTRA_CONFIGS``. A build with ``CONFIGS`` or ``EXTRA_CONFIGS`` names its executables without a configuration (with ``NEW_BINARY_NAMING``): ``seissol``, ``proxyseissol`` and ``testseissol-<module>``, followed by ``-<CUSTOM_BINARY_SUFFIX>`` if that is set. A build of a single configuration names them after it, e.g. ``seissol-viscoelastic-m3-o6-f64``.
- ``NEW_BINARY_NAMING``: (default on) use an updated binary naming scheme which uses lower-case letters and a shorter notation

Besides these, the host or, if enabled, the device architecture and backend are also encoded in the name of the executable.

.. _several_configurations:

Several configurations in one run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(NOTE: this feature is still considered highly experimental; its interface may change without notice)

The cells of a mesh group can compute in another configuration built into the executable than the rest of the mesh. ``ConfigMap`` in the ``equations`` section assigns them, e.g.

.. code-block:: Fortran

  &equations
  Configuration = 'elastic-o4-f64'
  ConfigMap = '2,3:viscoelastic-m3-o4-f64; 4:elastic-o5-f64'
  /

computes the groups 2 and 3 viscoelastically, the group 4 with order 5, and all other groups in the configuration of the run.

- The material file is queried for the cells of each configuration separately, with the material of that configuration. A group only has to define the parameters of the material it computes with.
- At a face between cells of different configurations, each cell reads the time integral of its neighbor converted into its own configuration: quantities with the same name are taken over, all others are zero; the precision is converted; a neighbor of a lower order is padded with zeros, and the trace of a neighbor of a higher order is projected to the face basis of the cell. So both cells see the exact traces of each other. The neighbor enters the Riemann problem with the parameters it is posed with. A solid and a fluid pose the Riemann problem in different quantities; there, the conversion gives the part of the neighbor in the state on the face, in the quantities of the cell, with a matrix per face that the initialization sets up from the impedances of both sides.
- Such faces are possible between configurations whose materials pose the Riemann problem in the same quantities (elastic, viscoelastic and anisotropic; acoustic and viscoacoustic; poroelastic only with itself), between an isotropic solid (elastic, viscoelastic) and a fluid (acoustic, viscoacoustic), and that fuse the same number of simulations. Between an anisotropic and an isotropic solid, the Riemann problem is posed in the anisotropic material, into which the isotropic one converts exactly (a viscoelastic one with the elastic part it poses its Riemann problem with). At a face between a solid and a fluid, both meet in the traction on the face and the velocity, with the impedance of each side: the fluid has none against shear, so that the shear traction on the face vanishes, and the pressure of the fluid and the normal stress of the solid meet as the normal traction. That is the same face as one between a solid and a fluid with zero shear modulus within the elastic equations. All configurations of a run have to fuse the same number of simulations. On GPUs, such faces need TensorForge as the device code generator (``DEVICE_CODEGEN``). Dynamic rupture faces between configurations and checkpoints of such runs are not supported yet.
- A cell weighs in the partitioning and the LTS clustering with the cost of its configuration relative to the configuration of the run, as Mini SeisSol measures it (see :ref:`env_vars`).
- The initial condition is set up for each configuration. Configurations with the same material set it up alike, from the material of the same cell. A planar wave is therefore an exact solution across configurations of the same material only.
- The wave field and the free surface are written into one output for all configurations. A quantity is written once for the cells of all configurations that have it, and as NaN for the others; the cell data ``config`` gives the position of the configuration of a cell in the list of the run: 0 for ``Configuration``, then the ones of ``ConfigMap`` in the order of the build. The log lists the configurations of the run. The analysis output writes one file per configuration.

Generic parameters
------------------

- ``LOG_LEVEL_MASTER``: determines the minimum relevance level of log messages that are printed. Possible values are ``error``, ``warning``, ``info``, and ``debug``.
- ``ASAGI``: enable or disable ASAGI as input for fault and material information. Note that easi will need to be built with ASAGI support if you select this option.
- ``COVERAGE``: determine the code coverage (currently only relevant for the CI).
- ``DR_QUAD_RULE``: the quadrature rule used for the Dynamic Rupture, the same for all configurations of the executable; the help of the executable states it (``seissol --help``). Currently permits either ``dunavant`` (the default) or ``stroud``, referring to the respective quadrature scheme names. ``dunavant`` uses Dunavant's symmetric rules exact up to degree 2 · order, with 6, 12, 16, 25, 33, 42 and 61 points for the orders 2 to 8 (degree 17 for order 8, since the rule of degree 16 has points outside the triangle); ``stroud`` is the conical product rule with (order + 1)² points, exact up to degree 2 · order + 1, whose points crowd towards one corner of the face. Both integrate the product of two basis functions of a face exactly; with its fewer points, ``dunavant`` takes about 30 % less work per fault face at order 6.
- ``PLASTICITY_METHOD``: the points at which the plastic correction is computed, and with them the plasticity matrices to be used. Options are

    * ``sq``: symmetric quadrature (the default), the points of a quadrature rule that integrates the product of two basis functions exactly, is invariant under all permutations of the vertices and has all points inside the cell (4, 14, 24, 46, 79, 128 and 188 points for the orders 2 to 8). Converting back to the basis is the L2 projection with that rule.
    * ``nb``: nodal basis, the warp & blend nodes of the basis: one point per basis function, some of them on the boundary of the cell. Converting back to the basis interpolates.
    * ``ip``: interpolation points, the conical-product (Stroud) quadrature rule with (order + 1)^3 points. Converting back to the basis is the L2 projection with that rule.
- ``PROXY_PYBINDING``: compile Python bindings for the SeisSol proxy
- ``TESTING``: compile unit tests for SeisSol
- ``TESTING_GENERATED``: compile unit tests for code-generated kernels. Useful for either developing a code generator, or when figuring out why a simulation runs into Inf/NaN or invalid results. Without it, the code generator writes no unit tests either.
- ``TESTING_PYTHON``: run additional code generator Python unit tests (those can also be run separately using ``pytest`` directly)
- ``SEISSOL_VERIFICATION_CONFIGS``: the configurations of the executable that the verification cases of the tests (``tests/verification``, end-to-end runs on generated mini meshes) run in: ``all`` of them (default), ``first``, or a list of their names. A case of a configuration is named after it, e.g. ``verification.smoke-box48<elastic-o6-f64>``.
- ``SANITIZERS``: instrument the C and C++ code of the build (SeisSol, its submodules, the generated kernels and the tests; not the device code) with the given sanitizers, as a list of ``address``, ``undefined``, ``leak``, ``thread`` and ``memory``, e.g. ``address;undefined``. ``thread`` needs Clang with the OpenMP runtime of LLVM, ``memory`` Clang and dependencies built with it as well. The tests fail on any report of a sanitizer; ``.ci/sanitizers/environment.sh`` sets the runtime options the CI runs them with. Combine it with ``CMAKE_BUILD_TYPE=Debug`` and ``-O1 -g1`` as ``CMAKE_C_FLAGS_DEBUG`` and ``CMAKE_CXX_FLAGS_DEBUG``, which keeps the assertions and the line numbers of the reports (with full debug information, GCC compiles the generated kernels many times as long), and with ``GEMM_TOOLS_LIST=none``, since no sanitizer sees into the kernels of PSpaMM or LIBXSMM.
- ``LTO``: enable link-time optimization
- ``GRAPH_PARTITIONING_LIBS``: compile for the given graph partitioning libraries. Allowed options are:

    * ``parmetis``: repository: https://github.com/KarypisLab/ParMETIS
    * ``PTScotch``: repository: https://www.labri.fr/perso/pelegrin/scotch/
    * ``ParHIP``: repository: https://github.com/KaHIP/KaHIP
    * ``none``: do not require a graph partitioning library. Recommended only for single-node/proxy builds.
- ``BUILD_DOCS``: build the (sphinx) documentation as HTML
- ``BUILD_DOXYGEN``: build the doxygen docs for the SeisSol code

CPU-specific parameters
-----------------------

- ``HOST_ARCH``: the parameter to tune the architecture for. See build-archs for an overview. Applies both compiler options and GEMM generation. If you get an illegal instruction error, check this variable first. Since there are many options, we have an own page for them. In short, use ``hsw`` when using a Linux/Windows PC or an older Mac. Use ``apple-m1`` or ``apple-m2`` on the latest Mac computers.
- ``MEMKIND``: enables HBM support.
- ``NUMA_AWARE_PINNING``: pin the free CPUs (those used for the communication and IO threads) according to the given NUMA domains.
- ``MEMORY_LAYOUT``: the sparsity patterns to apply. If set to ``auto``, code generation chooses
  a layout from the target-specific config folder (``codegen/config/cpu`` or
  ``codegen/config/gpu``). Passing a filename such as ``dense.xml`` or ``ms32.xml`` resolves it
  from that same target-specific folder automatically.
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

    * ``none``: compile for the CPU
    * ``cuda``: Nvidia CUDA
    * ``hip``: AMD HIP, using ROCm (or CUDA)
    * ``acpp``: SYCL, more specifically AdaptiveCpp, formerly known as Open SYCL and hipSYCL. Provides support for Intel, AMD, and Nvidia GPUs. Repository: https://github.com/AdaptiveCpp/AdaptiveCpp
    * ``oneapi``: SYCL, more specifically Intel Data Parallel C++ (DPC++). Provides support for Intel, AMD, and Nvidia GPUs. The open source variant is located under https://github.com/intel/llvm
- ``DEVICE_ARCH``: the parameter to tune and compile the kernels for. See build-archs for an overview.
- ``DEVICE_CODEGEN``: the GPU code generator to use. Unlike the CPU code generator list, you may only choose a single code generator here.

    * ``gemmforge-chainforge``: Gemmforge and chainforge (auto-enabled if the latter is found)
    * ``tensorforge``: TensorForge
    * ``tinytc``: Tiny Tensor Compiler (for Intel)
- ``DEVICE_CODEGEN_SHARDS``: the number of files the generated GPU routines are spread over, so that they compile at the same time; ``0`` (the default) takes four per configuration. The routines are spread by size, so the files take about the same time.
- ``DEVICE_KERNEL_INFOPRINT``: print register/resource usage info for each compiled kernel while compiling SeisSol.
- ``DEVICE_KERNEL_SAVETEMPS``: store the temporary output (e.g. PTX, AMDGCN ISA) of the GPU compilers. Useful e.g. for kernel debugging.
- ``SYCL_USE_NVHPC``: if AdaptiveCpp is compiled with NVHPC support, and we use NVHPC
- ``USE_GRAPH_CAPTURING``: if a compute graph feature is available, then use it. This is currently the case for CUDA (since 11.0), HIP (requires ROCm 6.1 or higher), or oneAPI.
- ``ENABLE_PROFILING_MARKERS``: Currently available for CUDA and HIP
- ``DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS``: enable some faster, handwritten accumulation kernels (default for NVIDIA and AMD)
