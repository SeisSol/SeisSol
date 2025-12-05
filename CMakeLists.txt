# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

cmake_minimum_required(VERSION 3.20)

# use <PackageName>_ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
# honor CMAKE_REQUIRED_LIBRARIES in the include file check macros
if(POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

project(SeisSol LANGUAGES C CXX)

# TODO: Move one dir up.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake"
                      "${PROJECT_SOURCE_DIR}/submodules/Device/cmake"
                      ${CMAKE_MODULE_PATH})

if (CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "22.1.0")
    message(FATAL_ERROR "NVHPC version 22.1.0 or higher is required")
  endif()
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "24.5.0")
    message(WARNING "NVHPC version 24.7.0 or higher is recommended")
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(GNU_REQUIRED_VERSION "9.0")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS GNU_REQUIRED_VERSION)
    message(WARNING "GCC versions below ${GNU_REQUIRED_VERSION} are not supported.")
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  message(WARNING "The Intel Compiler Classic (ICC) is deprecated; consider switching to the Intel oneAPI C++ Compiler (ICX) instead.")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Cray|CrayClang")
  message(WARNING "The communication thread is not supported for the Cray compiler. Make sure to set `SEISSOL_COMMTHREAD=0` when running SeisSol in multi-node setups.")
endif()

# find Python early
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# set hardware specific definition needed for seissol compilation
# 'process_users_input' returns the following:
#
#       switches: HDF5, NETCDF, GRAPH_PARTITIONING_LIBS, MPI, OPENMP, ASAGI, MEMKIND,
#                 PROXY_PYBINDING, ENABLE_PIC_COMPILATION, PREMULTIPLY_FLUX, SHARED
#
#       user input:   HOST_ARCH, DEVICE_ARCH, DEVICE_SUB_ARCH,
#                     ORDER, NUMBER_OF_MECHANISMS, EQUATIONS,
#                     PRECISION, DYNAMIC_RUPTURE_METHOD,
#                     NUMBER_OF_FUSED_SIMULATIONS,
#                     MEMORY_LAYOUT,
#                     LOG_LEVEL, LOG_LEVEL_MASTER,
#                     GEMM_TOOLS_LIST, EXTRA_CXX_FLAGS,
#                     DR_QUAD_RULE
#       derived input: REAL_SIZE_IN_BYTES, ALIGNMENT, VECTORSIZE ARCH_STRING, WITH_GPU, DEVICE_VENDOR
#
include(cmake/process_users_input.cmake)

if (TESTING AND COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
  message(STATUS "Enable coverage computation.")
  include(cmake/CodeCoverage.cmake)
  append_coverage_compiler_flags()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
  message(STATUS "Set build type to Release as none was supplied.")
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Release" "Debug" "RelWithDebInfo") # MinSizeRel is useless for us

if(PROXY_PYBINDING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

find_package(GemmTools REQUIRED)

add_subdirectory(codegen)

add_library(seissol-common-properties INTERFACE)

if (WITH_GPU)
  set(SEISSOL_DEVICE_INCLUDE ${DEVICE_INCLUDE_DIRS})

  if (DEVICE_BACKEND STREQUAL "cuda")
    include(cmake/cuda.cmake)
  endif()
  if (DEVICE_BACKEND STREQUAL "hip")
    include(cmake/hip.cmake)
  endif()
  if (DEVICE_BACKEND STREQUAL "hipsycl" OR DEVICE_BACKEND STREQUAL "acpp" OR DEVICE_BACKEND STREQUAL "oneapi")
    include(cmake/sycl.cmake)
  endif()
endif()

target_compile_options(seissol-common-properties INTERFACE ${EXTRA_CXX_FLAGS})
target_compile_options(seissol-common-properties INTERFACE -fPIC)

add_subdirectory(src)

# enforce code generation to run before any other target
add_dependencies(seissol-common-properties seissol-codegen)

target_link_libraries(seissol-lib PUBLIC seissol-common-properties)
target_link_libraries(seissol-kernel-lib PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PRIVATE seissol-kernel-lib)

if (WITH_GPU)
  add_dependencies(seissol-device-lib seissol-codegen)
endif()

if(GemmTools_INCLUDE_DIRS)
  target_include_directories(seissol-common-properties INTERFACE ${GemmTools_INCLUDE_DIRS})
endif()
if(GemmTools_LIBRARIES)
  target_link_libraries(seissol-common-properties INTERFACE ${GemmTools_LIBRARIES})
endif()
if(GemmTools_COMPILER_DEFINITIONS)
  target_compile_definitions(seissol-common-properties INTERFACE ${GemmTools_COMPILER_DEFINITIONS})
endif()

# Find appropriate compiler flags based on the target computer architecture
# and the target compiler
# 'get_arch_specific_flags' returns: CPU_ARCH_FLAGS
include(cmake/cpu_arch_flags.cmake)
get_arch_flags(${HOST_ARCH} ${CMAKE_CXX_COMPILER_ID})

# set hardware/compiler specific definitions and flags
target_compile_options(seissol-common-properties INTERFACE ${CPU_ARCH_FLAGS})

target_compile_definitions(seissol-common-properties INTERFACE LOGLEVEL=${LOG_LEVEL})
target_compile_definitions(seissol-common-properties INTERFACE LOG_LEVEL=${LOG_LEVEL_MASTER}
                                                     LOGLEVEL0=${LOG_LEVEL_MASTER})

# enable interproc. opts for small cores
#if cpu in ['knc', 'knl', 'skx']:
#    flags.extend(['-ip'])

# Libs
include(ExternalProject)

find_package(easi 1.0.0 REQUIRED)
target_link_libraries(seissol-common-properties INTERFACE easi::easi)

if (OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(seissol-common-properties INTERFACE OpenMP::OpenMP_CXX)
  target_compile_definitions(seissol-common-properties INTERFACE OMP OMPI_SKIP_MPICXX)
endif()

if (MPI)
  find_package(MPI REQUIRED)

  target_include_directories(seissol-common-properties SYSTEM INTERFACE ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(seissol-common-properties INTERFACE MPI::MPI_C)

  target_compile_definitions(seissol-common-properties INTERFACE USE_MPI PARALLEL)
endif()

if (NUMA_AWARE_PINNING)
  target_compile_definitions(seissol-common-properties INTERFACE USE_NUMA_AWARE_PINNING)
  find_package(NUMA REQUIRED)

  target_link_libraries(seissol-common-properties INTERFACE NUMA::NUMA)
endif()

if (ASAGI)
  find_package (asagi)
  if (NOT asagi_FOUND)
      find_package (PkgConfig REQUIRED)
      pkg_check_modules (ASAGI REQUIRED IMPORTED_TARGET asagi)
      add_library(asagi::asagi ALIAS PkgConfig::ASAGI)
  endif()

  target_link_libraries(seissol-common-properties INTERFACE asagi::asagi)
  target_compile_definitions(seissol-common-properties INTERFACE USE_ASAGI)
endif()

#set(HDF5_PREFER_PARALLEL True)
if (NETCDF)
  find_package(NetCDF REQUIRED)
  target_include_directories(seissol-common-properties INTERFACE ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(seissol-common-properties INTERFACE ${NetCDF_LIBRARY})
  target_compile_definitions(seissol-common-properties INTERFACE USE_NETCDF)
endif()

if (HDF5)
  if (MPI)
    set(HDF5_PREFER_PARALLEL True)
  endif()
  if (NETCDF)
    find_package(HDF5 REQUIRED COMPONENTS C HL)
    target_include_directories(seissol-common-properties INTERFACE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(seissol-common-properties INTERFACE ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})
  else()
    find_package(HDF5 REQUIRED COMPONENTS C)
    target_include_directories(seissol-common-properties INTERFACE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(seissol-common-properties INTERFACE ${HDF5_C_LIBRARIES})
  endif()
  target_compile_definitions(seissol-common-properties INTERFACE USE_HDF)
  if (NOT HDF5_IS_PARALLEL)
    message(WARNING "The found parallel Hdf5 installation is not parallel. Compile errors will occur.")
  endif()
endif()

# Parmetis
if ("parmetis" IN_LIST GRAPH_PARTITIONING_LIBS)
  find_package(METIS REQUIRED)
  if (NOT METIS_64_BIT_INTEGER)
      message(WARNING "Found METIS compiled with IDXTYPEWIDTH = 32. It is strongly recommend to compile METIS with IDXTYPEWIDTH = 64, because otherwise the partitioning of large meshes might fail.")
  endif()
  find_package(ParMETIS REQUIRED)
  target_link_libraries(seissol-common-properties INTERFACE ParMETIS::ParMETIS)
  target_compile_definitions(seissol-common-properties INTERFACE USE_METIS USE_PARMETIS)
endif()
if ("parhip" IN_LIST GRAPH_PARTITIONING_LIBS)
  find_package(ParHIP REQUIRED)
  target_link_libraries(seissol-common-properties INTERFACE ParHIP::ParHIP)
  target_compile_definitions(seissol-common-properties INTERFACE USE_PARHIP)
endif()
if ("ptscotch" IN_LIST GRAPH_PARTITIONING_LIBS)
  find_package(SCOTCH REQUIRED)
  target_link_libraries(seissol-common-properties INTERFACE SCOTCH::ptscotch SCOTCH::ptscotcherr)
  target_compile_definitions(seissol-common-properties INTERFACE USE_PTSCOTCH)
endif()
if (NOT GRAPH_PARTITIONING_LIBS)
  message(WARNING "Compiling without graph partitioning library; expect poor performance in multi-rank runs.")
endif()

if (MEMKIND)
  find_package(Memkind REQUIRED)
  target_include_directories(seissol-common-properties INTERFACE ${MEMKIND_INCLUDE_DIR})
  target_link_libraries(seissol-common-properties INTERFACE ${MEMKIND_LIBRARIES})
  target_compile_definitions(seissol-common-properties INTERFACE USE_MEMKIND)
endif()

if (INTEGRATE_QUANTITIES)
  target_compile_definitions(seissol-common-properties INTERFACE INTEGRATE_QUANTITIES)
endif()

if (USE_DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
  target_compile_definitions(seissol-common-properties INTERFACE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
endif()

if (ADDRESS_SANITIZER_DEBUG)
  target_link_libraries(seissol-common-properties INTERFACE debug
          -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address -static-libasan
  )
endif()

find_package(Eigen3 REQUIRED)
target_link_libraries(seissol-common-properties INTERFACE Eigen3::Eigen)

# needed explicitly also to silence some clang-tidy differences
if (NOT WITH_GPU)
  # but: it still causes issues with some GPU builds until we refactor the DR parameter datastructures
  target_compile_definitions(seissol-common-properties INTERFACE EIGEN_INITIALIZE_MATRICES_BY_ZERO)
endif()

# the following check is needed because allowing 3.4 and 5.0.0 at the same
# time was not possible inline in the find_package
if ("${Eigen3_VERSION}" VERSION_LESS "3.4")
  message(FATAL_ERROR "Your installed Eigen version is too old. Needed: at least 3.4 .")
endif()

# needed because Eigen (also 5.0) still incorrectly computes some values and fails some unit tests
if(ENABLE_MATERIAL_POROELASTIC)
  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    include(FortranCInterface)
    FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
    find_package(LAPACK REQUIRED)
    target_include_directories(seissol-common-properties INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(seissol-common-properties INTERFACE ${LAPACK_LIBRARIES})
    target_compile_definitions(seissol-common-properties INTERFACE USE_LAPACK)
  else()
    message(FATAL_ERROR "SeisSol needs a Fortran compiler.")
  endif()
endif()

find_package(yaml-cpp REQUIRED)

if ("${yaml-cpp_VERSION}" VERSION_GREATER_EQUAL "0.8")
  target_link_libraries(seissol-common-properties INTERFACE yaml-cpp::yaml-cpp)
elseif ("${yaml-cpp_VERSION}" VERSION_GREATER_EQUAL "0.7")
  target_link_libraries(seissol-common-properties INTERFACE yaml-cpp)
else()
  # fallback code for old versions
  if (YAML_CPP_INCLUDE_DIR AND EXISTS "${YAML_CPP_INCLUDE_DIR}")
    target_include_directories(seissol-common-properties INTERFACE ${YAML_CPP_INCLUDE_DIR})
  endif()
  if (YAML_CPP_LIBRARIES)
    # use the YAML_CPP_LIBRARIES, if available (though it may just say `yaml-cpp`)
    target_link_libraries(seissol-common-properties INTERFACE ${YAML_CPP_LIBRARIES})
  else()
    # fallback
    target_link_libraries(seissol-common-properties INTERFACE yaml-cpp)
  endif()
endif()

# Note: eigen3 enables cuda support in its headers by default.
# The problem happens in `eigen3/Eigen/src/Core/util/Meta.h` while
# including `math_constants.h`. This header file is located in
# a different sub-folder in case of ROCm HIP. SeisSol does not
# use eigen in GPU kernels. Therefore, it is better to change
# the default behaviour (which assumes CUDA support) and
# set `EIGEN_NO_CUDA`
target_compile_definitions(seissol-common-properties INTERFACE EIGEN_NO_CUDA)

if (EIGEN3_INCLUDE_DIR)
  target_include_directories(seissol-common-properties INTERFACE ${EIGEN3_INCLUDE_DIR})
endif()

target_include_directories(seissol-common-properties INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src/
    ${CMAKE_BINARY_DIR}/codegen/
)

# submodule includes

# Note: it is better to include `async` as
# system headers because they emit lots of warnings
# from clang. Most of them are issues with respect
# to overriden virtual methods

# for the same reason, do the same for all other submodules (for now)

# not explicitly listed: PUML, utils, fty, xdmfwriter
target_include_directories(seissol-common-properties SYSTEM INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/submodules
    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/async
)

if (PREMULTIPLY_FLUX)
  target_compile_definitions(seissol-common-properties INTERFACE USE_PREMULTIPLY_FLUX)
endif()

if (NOT NEW_BINARY_NAMING)
  # adjust prefix name of executables
  if ("${DEVICE_ARCH}" STREQUAL "none")
    set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${PRECISION_PREFIX}${HOST_ARCH}_${ORDER}_${EQUATIONS}")
  else()
    set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${PRECISION_PREFIX}${DEVICE_ARCH}_${DEVICE_BACKEND}_${ORDER}_${EQUATIONS}")
  endif()

  if(${NUMBER_OF_FUSED_SIMULATIONS} GREATER 1)
    set(EXE_NAME_PREFIX "${EXE_NAME_PREFIX}_f${NUMBER_OF_FUSED_SIMULATIONS}")
  endif()
else()
  if (WITH_GPU)
    set(NAME_SUFFIX "gpu")
  else()
    set(NAME_SUFFIX "cpu")
  endif()
  if (NOT MULTISETUP)
    if (EQUATIONS STREQUAL "viscoelastic2")
      set(NAME_SUFFIX "${NAME_SUFFIX}-viscoelastic")
    else()
      set(NAME_SUFFIX "${NAME_SUFFIX}-${EQUATIONS}")
    endif()
    if (NUMBER_OF_MECHANISMS GREATER 0)
      set(NAME_SUFFIX "${NAME_SUFFIX}-${NUMBER_OF_MECHANISMS}")
    endif()
    set(NAME_SUFFIX "${NAME_SUFFIX}-p${ORDER}")
    if (PRECISION STREQUAL "single")
      set(NAME_SUFFIX "${NAME_SUFFIX}-f32")
    else()
      set(NAME_SUFFIX "${NAME_SUFFIX}-f64")
    endif()
    if(${NUMBER_OF_FUSED_SIMULATIONS} GREATER 1)
      set(NAME_SUFFIX "${NAME_SUFFIX}-s${NUMBER_OF_FUSED_SIMULATIONS}")
    endif()
  endif()
endif()

if (SHARED)
  if (NOT NEW_BINARY_NAMING)
    set_target_properties(seissol-lib PROPERTIES OUTPUT_NAME "SeisSol_lib_${EXE_NAME_PREFIX}")
  else()
    set_target_properties(seissol-lib PROPERTIES OUTPUT_NAME "seissollib-${NAME_SUFFIX}")
  endif()
  install(TARGETS seissol-lib LIBRARY)
endif()

if (WITH_GPU)
  if (DEVICE_BACKEND STREQUAL "cuda")
    target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_CUDA)
  endif()
  if (DEVICE_BACKEND STREQUAL "hip")
    target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_HIP)
  endif()
  if (DEVICE_BACKEND STREQUAL "hipsycl" OR DEVICE_BACKEND STREQUAL "acpp" OR DEVICE_BACKEND STREQUAL "oneapi")
    target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_SYCL)
  endif()

  # set SeisSol GPU definitions
  target_compile_definitions(seissol-common-properties INTERFACE ACL_DEVICE)

  # add device submodule
  add_subdirectory(submodules/Device)
  target_link_libraries(seissol-common-properties INTERFACE device)

  # add SeisSol GPU part
  target_link_libraries(seissol-common-properties INTERFACE seissol-device-lib)

  if (NOT NEW_BINARY_NAMING)
    set_target_properties(seissol-device-lib PROPERTIES OUTPUT_NAME "SeisSol_gpucodegen_${EXE_NAME_PREFIX}")
    set_target_properties(device PROPERTIES OUTPUT_NAME "device_${CMAKE_BUILD_TYPE}_${DEVICE_BACKEND}")
  else()
    set_target_properties(seissol-device-lib PROPERTIES OUTPUT_NAME "seissol-gpucode-${NAME_SUFFIX}")
  endif()

  # set up GPU install targets
  install(TARGETS seissol-device-lib LIBRARY)
  install(TARGETS device LIBRARY)
endif()

# C++ compiler settings
# TODO: mostly remove or make optional? (except the EIGEN_DONT_VECTORIZE for NVHPC below)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(seissol-common-properties INTERFACE $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-pedantic -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas>)

  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  # todo remove std?, is ffreestanding needed?
  target_compile_options(seissol-common-properties INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-pedantic -std=c++17>)

  # Activate interprocedual optimization.
  #set_property(TARGET seissol-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
  target_compile_options(seissol-common-properties INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-pedantic -Wall -Wextra -Wno-unused-parameter>)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
    # warning 186 - pointless comparison of unsigned integer with zero detected
    # warning 612 - overloaded virtual function is only partially overridden
    # warning 111 - statement is unreachable
    target_compile_options(seissol-common-properties INTERFACE $<$<COMPILE_LANGUAGE:CXX>:--display_error_number --diag_suppress186 --diag_suppress612 --diag_suppress111>)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "22.3.1")
      target_compile_options(seissol-common-properties INTERFACE $<$<COMPILE_LANGUAGE:CXX>:--diag_suppress1>)
    endif()

    # NOTE:
    # 1. --pending_instantiations=0 allows an infinite recursive template instantiation
    target_compile_options(seissol-common-properties INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wc,--pending_instantiations=0>)

    if ("${Eigen3_VERSION}" VERSION_LESS "5.0")
      # 2. EIGEN_DONT_VECTORIZE=1 ; Eigen3 3.4 does not have support for NVHPC vectorization
      target_compile_definitions(seissol-common-properties INTERFACE EIGEN_DONT_VECTORIZE=0)
    endif()
endif()

find_package(FILESYSTEM REQUIRED)
target_link_libraries(seissol-common-properties INTERFACE std::filesystem)

# put every executable and shared library in the root build directory at the end
# (for compatibility reasons; TODO: remove at some point)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (LIKWID)
  find_package(likwid REQUIRED)
  target_compile_definitions(seissol-lib PUBLIC LIKWID_PERFMON)
  target_include_directories(seissol-lib PUBLIC ${LIKWID_INCLUDE_DIR})
  target_link_libraries(seissol-lib PUBLIC likwid::likwid)
endif()

add_subdirectory(external)

if (TESTING)
  enable_testing()
  include(cmake/doctest.cmake)

  add_subdirectory(tests)

  # Coverage
  if(COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
    setup_target_for_coverage_lcov(
            NAME seissol-coverage
            EXECUTABLE seissol-serial-test
            EXCLUDE "/usr/*"
                    "submodules/*"
                    "tests/*"
                    "external/*"
                    "*/codegen/GeneratedCode/*"
                    "/opt/*"
            LCOV_ARGS --ignore-errors inconsistent
                      --ignore-errors unused
                      --ignore-errors mismatch
    )
  endif()
endif()

add_subdirectory(app)

if (WITH_GPU)
  string(TOUPPER "${DEVICE_BACKEND}" BACKEND_UPPER_CASE)
  set(RUN_GPU_TIME_SELECTOR "${BACKEND_UPPER_CASE}_VISIBLE_DEVICES")
  configure_file(${PROJECT_SOURCE_DIR}/cmake/gpu-launch-template.in
                 ${PROJECT_BINARY_DIR}/seissol-launch
                 @ONLY)
  file(CHMOD ${PROJECT_BINARY_DIR}/seissol-launch PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
  install(FILES ${PROJECT_BINARY_DIR}/seissol-launch DESTINATION ${CMAKE_INSTALL_PREFIX}/share
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
endif()

if (BUILD_DOXYGEN)
  add_subdirectory(docs/doxygen)
endif()

if (BUILD_DOCS)
  add_subdirectory(docs)
endif()

# https://blog.kitware.com/static-checks-with-cmake-cdash-iwyu-clang-tidy-lwyu-cpplint-and-cppcheck/
# https://ortogonal.github.io/cmake-clang-tidy/
