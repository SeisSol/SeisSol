# SPDX-FileCopyrightText: 2019-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

cmake_minimum_required(VERSION 3.13)

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
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "24.7.0")
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

set(HARDWARE_DEFINITIONS "ALIGNMENT=${ALIGNMENT}"
                         "VECTORSIZE=${VECTORSIZE}"
                         "REAL_SIZE=${REAL_SIZE_IN_BYTES}"
                         "ALIGNED_REAL_SIZE=${REAL_SIZE_IN_BYTES}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
  message(STATUS "Set build type to Release as none was supplied.")
endif()

if (NOT (DEVICE_BACKEND STREQUAL "hip" AND DEVICE_ARCH MATCHES "sm_*"))
  # Note: hip 4.x has a bug during setting rpath while compiling
  #       device code for nvidia gpus. This scenario is used in
  #       SeisSol's GPU CI pipeline.
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
endif()

if(NUMBER_OF_FUSED_SIMULATIONS GREATER 1)
  message(FATAL_ERROR "SeisSol does not support Fused Simulations yet")
endif()

# LTO
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${LTO})

# Generate version.h
include(GetGitRevisionDescription)

# get GIT info
git_describe(PACKAGE_GIT_VERSION --tags --always --dirty=\ \(dirty\) --broken=\ \(broken\))
if (${PACKAGE_GIT_VERSION} MATCHES "NOTFOUND")
  set(PACKAGE_GIT_VERSION "(unknown)")
  set(PACKAGE_GIT_HASH "(unknown)")
  set(PACKAGE_GIT_TIMESTAMP "9999-12-31T00:00:00+00:00")
else()
  get_git_commit_info(PACKAGE_GIT_HASH PACKAGE_GIT_TIMESTAMP)
endif()
string(SUBSTRING ${PACKAGE_GIT_TIMESTAMP} 0 4 PACKAGE_GIT_YEAR)

# write file and print info
configure_file("src/VersionTemplate.h"
               "${CMAKE_CURRENT_BINARY_DIR}/src/Version.h")
message(STATUS "Version: " ${PACKAGE_GIT_VERSION})
message(STATUS "Last commit: ${PACKAGE_GIT_HASH} at ${PACKAGE_GIT_TIMESTAMP}")

add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/generated_code")

set(GENERATED_FILES_FOR_SEISSOL src/generated_code/subroutine.h
                                src/generated_code/tensor.cpp
                                src/generated_code/subroutine.cpp
                                src/generated_code/gpulike_subroutine.cpp
                                src/generated_code/tensor.h
                                src/generated_code/init.cpp
                                src/generated_code/init.h
                                src/generated_code/kernel.h
                                src/generated_code/kernel.cpp)

if (TESTING AND TESTING_GENERATED)
  set(GENERATED_FILES_FOR_SEISSOL ${GENERATED_FILES_FOR_SEISSOL}
                                src/generated_code/test-kernel.cpp)
endif()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/)

if(PROXY_PYBINDING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

if(PREMULTIPLY_FLUX)
  set(PREMULTIPLY_FLUX_VALUE "--enable_premultiply_flux")
else()
  set(PREMULTIPLY_FLUX_VALUE "--disable_premultiply_flux")
endif()

find_package(GemmTools REQUIRED)

add_custom_command(
  COMMAND
  "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/generated-code/generate.py"
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/generated-code/matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code
     "--host_arch" ${HOST_ARCH_STR}
     "--device_codegen" ${DEVICE_CODEGEN}
     "--device_arch" ${DEVICE_ARCH_STR}
     "--device_backend" ${DEVICE_BACKEND}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" ${MEMORY_LAYOUT}
     "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
     "--gemm_tools" ${GEMM_TOOLS_LIST}
     "--drQuadRule" ${DR_QUAD_RULE}
     "--executable_libxsmm=${Libxsmm_executable_PROGRAM}"
     "--executable_pspamm=${PSpaMM_PROGRAM}"
     ${PREMULTIPLY_FLUX_VALUE} # boolean flag
     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/generated-code
     DEPENDS
        build-time-make-directory
        generated-code/generate.py
        generated-code/kernels
        generated-code/matrices
        generated-code/config
    OUTPUT
        ${GENERATED_FILES_FOR_SEISSOL}
    COMMENT "Code generation for tensor operations..."
       )

add_custom_target(seissol-codegen ALL DEPENDS ${GENERATED_FILES_FOR_SEISSOL})
add_library(seissol-common-properties INTERFACE IMPORTED)

include(src/sources.cmake)

target_compile_options(seissol-common-properties INTERFACE -fPIC)

# enforce code generation to run before any other target
add_dependencies(seissol-common-properties seissol-codegen)

target_link_libraries(seissol-lib PUBLIC seissol-common-properties)
target_link_libraries(seissol-kernel-lib PUBLIC seissol-common-properties)
target_link_libraries(seissol-common-lib PUBLIC seissol-common-properties)

target_link_libraries(seissol-common-lib PRIVATE seissol-kernel-lib)
target_link_libraries(seissol-lib PRIVATE seissol-kernel-lib)
target_link_libraries(seissol-lib PRIVATE seissol-common-lib)

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
target_compile_definitions(seissol-common-properties INTERFACE ${HARDWARE_DEFINITIONS})
target_compile_options(seissol-common-properties INTERFACE ${CPU_ARCH_FLAGS})

target_compile_definitions(seissol-common-properties INTERFACE LOGLEVEL=${LOG_LEVEL})
target_compile_definitions(seissol-common-properties INTERFACE LOG_LEVEL=${LOG_LEVEL_MASTER}
                                                     LOGLEVEL0=${LOG_LEVEL_MASTER})


if (PLASTICITY_METHOD STREQUAL "ip")
  target_compile_definitions(seissol-common-properties INTERFACE USE_PLASTICITY_IP)
elseif (PLASTICITY_METHOD STREQUAL "nb")
  target_compile_definitions(seissol-common-properties INTERFACE USE_PLASTICITY_NB)
endif()


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

find_package(PkgConfig REQUIRED) 
if (ASAGI)
  # todo warn if netcdf is off
  pkg_check_modules(ASAGI REQUIRED IMPORTED_TARGET asagi) # asagi_nompi?
  target_compile_definitions(seissol-common-properties INTERFACE USE_ASAGI)
  target_link_libraries(seissol-common-properties INTERFACE ${ASAGI_LDFLAGS})
  target_include_directories(seissol-common-properties INTERFACE ${ASAGI_INCLUDE_DIRS})
  target_compile_options(seissol-common-properties INTERFACE ${ASAGI_CFLAGS} ${ASAGI_CFLAGS_OTHER})
endif()

if (MEMKIND)
  find_package(Memkind REQUIRED)
  target_include_directories(seissol-common-properties INTERFACE ${MEMKIND_INCLUDE_DIR})
  target_link_libraries(seissol-common-properties INTERFACE ${MEMKIND_LIBRARIES})
  target_compile_definitions(seissol-common-properties INTERFACE USE_MEMKIND)
endif()

if(${EQUATIONS} STREQUAL "poroelastic")
  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    include(FortranCInterface)
    FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
    find_package(LAPACK REQUIRED)
    target_include_directories(seissol-common-properties INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(seissol-common-properties INTERFACE ${LAPACK_LIBRARIES})
  else()
    message(FATAL_ERROR "SeisSol needs a Fortran compiler.")
  endif()
endif()


if (INTEGRATE_QUANTITIES)
  target_compile_definitions(seissol-common-properties INTERFACE INTEGRATE_QUANTITIES)
endif()

if (DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
  target_compile_definitions(seissol-common-properties INTERFACE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
endif()

if (ADDRESS_SANITIZER_DEBUG)
  target_link_libraries(seissol-common-properties INTERFACE debug
          -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address -static-libasan
  )
endif()

# Note: it is better to include `async` as
# system headers because they emit lot's of warnings
# from clang. Most of them are issues with respect
# to overriden virtual methods
target_include_directories(seissol-common-properties SYSTEM INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/submodules/async
)

find_package(Eigen3 3.4 REQUIRED)
target_link_libraries(seissol-common-properties INTERFACE Eigen3::Eigen)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/submodules
    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
    ${CMAKE_CURRENT_BINARY_DIR}/src/
)

target_compile_definitions(seissol-common-properties INTERFACE
    CONVERGENCE_ORDER=${ORDER}
    NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
    ${DR_QUAD_RULE}
)

if (PREMULTIPLY_FLUX)
  target_compile_definitions(seissol-common-properties INTERFACE USE_PREMULTIPLY_FLUX)
endif()

# adjust prefix name of executables
if ("${DEVICE_ARCH_STR}" STREQUAL "none")
  set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${HOST_ARCH_STR}_${ORDER}_${EQUATIONS}")
else()
  set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${DEVICE_ARCH_STR}_${DEVICE_BACKEND}_${ORDER}_${EQUATIONS}")
endif()

if (SHARED)
  set_target_properties(seissol-lib PROPERTIES OUTPUT_NAME "SeisSol_lib_${EXE_NAME_PREFIX}")
  install(TARGETS seissol-lib LIBRARY)
endif()

if (WITH_GPU)
  if (DEVICE_BACKEND STREQUAL "cuda")
    target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_CUDA)
  endif()
  if (DEVICE_BACKEND STREQUAL "hip")
    target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_HIP)
  endif()
  if (DEVICE_BACKEND STREQUAL "hipsycl" OR DEVICE_BACKEND STREQUAL "oneapi")
    target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_SYCL)
  endif()

  # set SeisSol GPU definitions
  target_compile_definitions(seissol-common-properties INTERFACE ACL_DEVICE)

  # add device submodule
  add_subdirectory(submodules/Device)
  target_include_directories(seissol-common-properties INTERFACE submodules/Device)
  target_link_libraries(seissol-common-properties INTERFACE device)

  # add SeisSol GPU part
  target_link_libraries(seissol-common-properties INTERFACE seissol-device-lib)

  set_target_properties(seissol-device-lib PROPERTIES OUTPUT_NAME "SeisSol_gpucodegen_${EXE_NAME_PREFIX}")
  set_target_properties(device PROPERTIES OUTPUT_NAME "device_${CMAKE_BUILD_TYPE}_${DEVICE_ARCH_STR}_${DEVICE_BACKEND}")

  # set up GPU install targets
  install(TARGETS seissol-device-lib LIBRARY)
  install(TARGETS device LIBRARY)
  target_sources(seissol-common-lib PRIVATE ${SYCL_ONLY_SRC_FILES})
endif()

# C++ compiler settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(seissol-lib INTERFACE -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas>)

  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  # todo remove std?, is ffreestanding needed?
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

  # Activate interprocedual optimization.
  #set_property(TARGET seissol-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE) 
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wno-unused-parameter")

elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
    # warning 186 - pointless comparison of unsigned integer with zero detected
    # warning 612 - overloaded virtual function is only partially overridden
    # warning 111 - statement is unreachable
    set(WARNINGS "--display_error_number --diag_suppress186 --diag_suppress612 --diag_suppress111")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "22.3.1")
      set(WARNINGS "${WARNINGS} --diag_suppress1")
    endif()

    # NOTE:
    # 1. --pending_instantiations=0 allows an infinite recursive template instantiation
    # 2. EIGEN_DONT_VECTORIZE=1 waiting for eigen3 support for nvhpc compiler collection w.r.t. vectorization
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wc,--pending_instantiations=0 ${WARNINGS} -DEIGEN_DONT_VECTORIZE=0")
endif()

find_package(FILESYSTEM REQUIRED)
target_link_libraries(seissol-common-properties INTERFACE std::filesystem)

# Generated code does only work without red-zone.
if (HAS_REDZONE)
  set_source_files_properties(
      ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone
  )
endif()

# put every executable in the root build directory at the end
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (LIKWID)
  find_package(likwid REQUIRED)
  target_compile_definitions(seissol-lib PUBLIC LIKWID_PERFMON)
  target_include_directories(seissol-lib PUBLIC ${LIKWID_INCLUDE_DIR})
  target_link_libraries(seissol-lib PUBLIC likwid::likwid)
endif()

if (TESTING)
  enable_testing()
  include(cmake/doctest.cmake)

  set(seissol_test_sources
          src/tests/Model/TestModel.cpp
          src/tests/Initializer/TestInitializer.cpp
          src/tests/Numerical/TestNumerical.cpp
          src/tests/Geometry/TestGeometry.cpp
          src/tests/Kernel/TestKernel.cpp
          src/tests/SourceTerm/TestSourceTerm.cpp
          src/tests/Parallel/TestParallel.cpp
          src/tests/ResultWriter/TestResultWriter.cpp
          src/tests/Solver/time_stepping/TestSolverTimeStepping.cpp
          src/tests/DynamicRupture/TestDynamicRupture.cpp
          src/tests/Common/TestCommon.cpp
          )


  if (TESTING_GENERATED)
    # FIXME: fix by removing the `using namespace` in the test files
    # and replace by ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/test-kernel.cpp
    file(GLOB_RECURSE KERNEL_TESTFILES ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/**/test-kernel.cpp)
    set(seissol_test_sources
            ${seissol_test_sources}
            ${KERNEL_TESTFILES}
            )
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Initializer/time_stepping/mesh.h5
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/mesh.h5
    COPYONLY
    )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Initializer/time_stepping/material.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/material.yaml
    COPYONLY
    )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/ResultWriter/receiver_correct.dat
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/receiver_correct.dat
    COPYONLY
  )

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/fsrm_source1.dat
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/fsrm_source1.dat
    COPYONLY
  )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/fsrm_source2.dat
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/fsrm_source2.dat
    COPYONLY
  )

  if (NETCDF)
    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/source_loh.nrf
            ${CMAKE_CURRENT_BINARY_DIR}/Testing/source_loh.nrf
            COPYONLY
    )
    set(seissol_test_sources ${seissol_test_sources} ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/ReaderTest.cpp)
  endif()

  add_executable(seissol-serial-test
            ${seissol_test_sources}
            src/tests/TestMain.cpp)
  target_link_libraries(seissol-serial-test PRIVATE seissol-lib)
  target_include_directories(seissol-serial-test PUBLIC external/)

  separate_arguments(TESTING_COMMAND_LIST NATIVE_COMMAND ${TESTING_COMMAND})
  set_property(TARGET seissol-serial-test PROPERTY CROSSCOMPILING_EMULATOR ${TESTING_COMMAND_LIST})
  doctest_discover_tests(seissol-serial-test)

  # Avoid duplicate definition of FLOP counters
  target_compile_definitions(seissol-serial-test PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)

  # Coverage
  if(COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
    setup_target_for_coverage_lcov(
            NAME seissol-coverage
            EXECUTABLE seissol-serial-test
            EXCLUDE "/usr/*"
                    "submodules/*"
                    "*/tests/*"
                    "external/*"
                    "*/src/generated_code/*"
                    "/opt/*"
            LCOV_ARGS --ignore-errors inconsistent
                      --ignore-errors unused
    )
  endif()
endif()

# build seissol-bin
add_executable(seissol-bin src/Main.cpp)
set_target_properties(seissol-bin PROPERTIES OUTPUT_NAME "SeisSol_${EXE_NAME_PREFIX}")

target_link_libraries(seissol-bin PUBLIC seissol-lib)

install(TARGETS seissol-bin RUNTIME)
# end build seissol-bin

# build seissol-proxy
add_subdirectory(auto-tuning/proxy/src)
# end build seissol-proxy

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

# https://blog.kitware.com/static-checks-with-cmake-cdash-iwyu-clang-tidy-lwyu-cpplint-and-cppcheck/
# https://ortogonal.github.io/cmake-clang-tidy/

