cmake_minimum_required(VERSION 3.10)

# use <PackageName>_ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
# honor CMAKE_REQUIRED_LIBRARIES in the include file check macros
if(POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

project(SeisSol LANGUAGES C CXX Fortran)

if (CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
  set(NVHPC_REQUIRED_VERSION "22.1.0")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS NVHPC_REQUIRED_VERSION)
    message(FATAL_ERROR "NVHPC version ${NVHPC_REQUIRED_VERSION} or higher is required")
  endif()
endif()

# set hardware specific definition needed for seissol compilation
# 'process_users_input' returns the following:
#
#       switches: HDF5, NETCDF, METIS, MPI, OPENMP, ASAGI, SIONLIB, MEMKIND, PROXY_PYBINDING, ENABLE_PIC_COMPILATION
#
#       user's input: HOST_ARCH, DEVICE_ARCH, DEVICE_SUB_ARCH,
#                     ORDER, NUMBER_OF_MECHANISMS, EQUATIONS,
#                     PRECISION, DYNAMIC_RUPTURE_METHOD,
#                     NUMBER_OF_FUSED_SIMULATIONS,
#                     MEMORY_LAYOUT, COMMTHREAD,
#                     LOG_LEVEL, LOG_LEVEL_MASTER,
#                     GEMM_TOOLS_LIST, EXTRA_CXX_FLAGS,
#                     DR_QUAD_RULE
#       derived input: REAL_SIZE_IN_BYTES, ALIGNMENT, ARCH_STRING, NUMBER_OF_QUANTITIES, WITH_GPU
#
include(cmake/process_users_input.cmake)
set(HARDWARE_DEFINITIONS "ALIGNMENT=${ALIGNMENT}"
                         "REAL_SIZE=${REAL_SIZE_IN_BYTES}"
                         "ALIGNED_REAL_SIZE=${REAL_SIZE_IN_BYTES}")

# TODO: Move one dir up.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${CMAKE_ROOT}/Modules")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()


# Generate version.h

include(GetGitRevisionDescription)
git_describe(PACKAGE_GIT_VERSION  --dirty=\ \(modified\) --always)
configure_file("src/version_template.h"
               "${CMAKE_CURRENT_BINARY_DIR}/src/version.h")
message(STATUS "Current Git description: " ${PACKAGE_GIT_VERSION})


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

add_custom_command(
  COMMAND
  "/usr/bin/env" python3 ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/generate.py
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code
     "--host_arch" ${HOST_ARCH_STR}
     "--device_arch" ${DEVICE_ARCH_STR}
     "--device_backend" ${DEVICE_BACKEND}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" ${MEMORY_LAYOUT}
     "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
     "--gemm_tools" ${GEMM_TOOLS_LIST}
     "--drQuadRule" ${DR_QUAD_RULE}
     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/generated_code
     DEPENDS
        build-time-make-directory
        generated_code/arch.py
        generated_code/memlayout.py
        generated_code/viscoelastic.py
        generated_code/DynamicRupture.py
        generated_code/Plasticity.py
        generated_code/multSim.py
        generated_code/aderdg.py
        generated_code/generate.py
        generated_code/Point.py
        generated_code/elastic.py
        generated_code/viscoelastic2.py
        generated_code/poroelastic.py
        generated_code/anisotropic.py
        generated_code/SurfaceDisplacement.py
        generated_code/NodalBoundaryConditions.py
        generated_code/matrices
    OUTPUT
        ${GENERATED_FILES_FOR_SEISSOL}
    COMMENT "Codegen for tensor stuff."
       )

add_custom_target(SeisSol-codegen ALL DEPENDS ${GENERATED_FILES_FOR_SEISSOL})

find_package(GemmTools REQUIRED)
include(src/sources.cmake)

# enforce code generation to run before any other target
add_dependencies(SeisSol-lib SeisSol-codegen)
if (WITH_GPU)
  add_dependencies(SeisSol-device-lib SeisSol-codegen)
endif()


if(GemmTools_INCLUDE_DIRS)
  target_include_directories(SeisSol-lib PUBLIC ${GemmTools_INCLUDE_DIRS})
endif()
if(GemmTools_LIBRARIES)
  target_link_libraries(SeisSol-lib PUBLIC ${GemmTools_LIBRARIES})
endif()
if(GemmTools_COMPILER_DEFINITIONS)
  target_compile_definitions(SeisSol-lib PUBLIC ${GemmTools_COMPILER_DEFINITIONS})
endif()

# Find appropriate compiler flags based on the target computer architecture
# and the target compiler
# 'get_arch_specific_flags' returns: CPU_ARCH_FLAGS
include(cmake/cpu_arch_flags.cmake)
get_arch_flags(${HOST_ARCH} ${CMAKE_CXX_COMPILER_ID})

# set hardware/compiler specific definitions and flags
target_compile_definitions(SeisSol-lib PUBLIC ${HARDWARE_DEFINITIONS})
target_compile_options(SeisSol-lib PUBLIC ${CPU_ARCH_FLAGS})

target_compile_definitions(SeisSol-lib PUBLIC LOGLEVEL=${LOG_LEVEL})
target_compile_definitions(SeisSol-lib PUBLIC LOG_LEVEL=${LOG_LEVEL_MASTER}
                                              LOGLEVEL0=${LOG_LEVEL_MASTER})


if (PLASTICITY_METHOD STREQUAL "ip")
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY_IP)
elseif (PLASTICITY_METHOD STREQUAL "nb")
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY_NB)
endif()


# enable interproc. opts for small cores
#if cpu in ['knc', 'knl', 'skx']:
#    flags.extend(['-ip'])

# Libs
include(ExternalProject)

find_package(easi 1.0.0 REQUIRED)
target_link_libraries(SeisSol-lib PUBLIC easi::easi)

if (OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(SeisSol-lib PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(SeisSol-lib PUBLIC OpenMP::OpenMP_Fortran)
  target_compile_definitions(SeisSol-lib PUBLIC OMP OMPI_SKIP_MPICXX)

  if (WITH_GPU)
    find_package(OpenMP_Offloading QUIET)
  else()
    # Add DR source files to SeisSol-lib if we compile without any GPU support
    target_sources(SeisSol-lib PUBLIC ${DYNAMIC_RUPTURE_SRC_FILES})
  endif()

  # Note: there may be a conflict between OpenMP offloading and libsxmm.
  # A collision can happen due to -march flags for CPUs and OpenMP GPU target
  # (e.g. with clang -Xopenmp-target -march=sm_75). Therefore, we separate
  # compilations for the source and generated codes and apply different
  # imported OpenMP targets to them.
  if (OpenMP_Offloading_FOUND)
    message(STATUS "Found OpenMP Offloading: TRUE")
    add_library(SeisSol-omp-offloaded SHARED ${DYNAMIC_RUPTURE_SRC_FILES} ${DYNAMIC_RUPTURE_OMP_OFFLOADING_SRC_FILES})
    target_link_libraries(SeisSol-omp-offloaded PRIVATE OpenMP::OpenMP_CXX_with_offloading)
    target_compile_definitions(SeisSol-omp-offloaded PRIVATE ACL_DEVICE_OFFLOAD)

    add_dependencies(SeisSol-omp-offloaded SeisSol-codegen)
  else()
    # Add DR source files to SeisSol-lib if the compiler doesn't support GPU offloading
    target_sources(SeisSol-lib PUBLIC ${DYNAMIC_RUPTURE_SRC_FILES})
  endif()
endif()

if (MPI)
  find_package(MPI REQUIRED)

  target_include_directories(SeisSol-lib SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC MPI::MPI_C)

  target_include_directories(SeisSol-lib SYSTEM PUBLIC ${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC MPI::MPI_Fortran)

  target_compile_definitions(SeisSol-lib PUBLIC USE_MPI PARALLEL)

  if (MINI_SEISSOL)
    target_compile_definitions(SeisSol-lib PUBLIC USE_MINI_SEISSOL)
  endif()
endif()

if (COMMTHREAD)
  target_compile_definitions(SeisSol-lib PUBLIC USE_COMM_THREAD)
endif()

if (NUMA_AWARE_PINNING)
  target_compile_definitions(SeisSol-lib PUBLIC USE_NUMA_AWARE_PINNING)
  find_package(NUMA REQUIRED)

  target_include_directories(SeisSol-lib SYSTEM PUBLIC ${NUMA_INCLUDE_DIR})
  target_link_libraries(SeisSol-lib PUBLIC ${NUMA_LIBRARY})
endif()

#set(HDF5_PREFER_PARALLEL True)
if (NETCDF)
  find_package(NetCDF REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC ${NetCDF_LIBRARY})
  target_compile_definitions(SeisSol-lib PUBLIC USE_NETCDF)
endif()

if (HDF5)
  if (MPI)
    set(HDF5_PREFER_PARALLEL True)
  endif()
  find_package(HDF5 REQUIRED
          COMPONENTS C HL)
  target_include_directories(SeisSol-lib PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_HDF)
endif()  

# Parmetis
if (METIS)
  find_package(METIS REQUIRED)
  if (NOT METIS_64_BIT_INTEGER)
    message(WARNING "Found METIS which does not support IDXTYPEWIDTH = 64. But this is required for SeisSol.")
  endif()
  find_package(ParMETIS REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${PARMETIS_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC ${PARMETIS_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_METIS)
endif()

find_package(PkgConfig REQUIRED) 
if (ASAGI)
  # todo warn if netcdf is off
  pkg_check_modules(ASAGI REQUIRED IMPORTED_TARGET asagi) # asagi_nompi?
  target_compile_definitions(SeisSol-lib PUBLIC USE_ASAGI)
  target_link_libraries(SeisSol-lib PUBLIC ${ASAGI_STATIC_LDFLAGS})
  target_include_directories(SeisSol-lib PUBLIC ${ASAGI_INCLUDE_DIRS})
  target_compile_options(SeisSol-lib PUBLIC ${ASAGI_CFLAGS} ${ASAGI_CFLAGS_OTHER})
endif()

if (MEMKIND)
  find_package(Memkind REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${MEMKIND_INCLUDE_DIR})
  target_link_libraries(SeisSol-lib PUBLIC ${MEMKIND_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_MEMKIND)
endif()

if(${EQUATIONS} STREQUAL "poroelastic")
  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    include(FortranCInterface)
    FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
    find_package(LAPACK REQUIRED)
    target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
    message(${LAPACK_LIBRARIES})
    target_link_libraries(SeisSol-lib PUBLIC ${LAPACK_LIBRARIES})
  else()
    message(FATAL_ERROR "SeisSol needs a Fortran compiler.")
  endif()
endif()


if (INTEGRATE_QUANTITIES)
  target_compile_definitions(SeisSol-lib PUBLIC INTEGRATE_QUANTITIES)
endif()

if (ADDRESS_SANITIZER_DEBUG)
  target_link_libraries(SeisSol-lib PUBLIC debug
          -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address
  )
endif()

find_package(yaml-cpp 0.6.2 REQUIRED)
target_include_directories(SeisSol-lib PUBLIC ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(SeisSol-lib PUBLIC ${YAML_CPP_LIBRARIES})

# Note: it is better to include `eigen3` `async` as
# system headers because they emit lot's of warnings
# from clang. Most of them are issues with respect
# to overriden virtual methods
target_include_directories(SeisSol-lib SYSTEM PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen3
  ${CMAKE_CURRENT_SOURCE_DIR}/submodules/async
)

target_include_directories(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/submodules
    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
    ${CMAKE_CURRENT_BINARY_DIR}/src/
)

target_compile_definitions(SeisSol-lib PUBLIC
    CONVERGENCE_ORDER=${ORDER}
    NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
    NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
    ENABLE_MATRIX_PREFETCH
    ${DR_QUAD_RULE}
)


if (WITH_GPU)
  add_subdirectory(submodules/Device)
  target_include_directories(SeisSol-lib PUBLIC submodules/Device)
  target_link_libraries(SeisSol-lib PUBLIC device)
  target_link_libraries(SeisSol-lib PUBLIC SeisSol-device-lib)
  target_compile_definitions(SeisSol-lib PUBLIC ACL_DEVICE)

  if (OpenMP_Offloading_FOUND)
    target_compile_definitions(SeisSol-lib PUBLIC ACL_DEVICE_OFFLOAD)

    get_target_property(SEISSOL_DEFENITIONS SeisSol-lib COMPILE_DEFINITIONS)
    target_compile_definitions(SeisSol-omp-offloaded PRIVATE ${SEISSOL_DEFENITIONS})

    get_target_property(SEISSOL_INCLUDES SeisSol-lib INCLUDE_DIRECTORIES)
    get_target_property(SEISSOL_SYSTEM_INCLUDES SeisSol-lib INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    list(REMOVE_ITEM SEISSOL_INCLUDES ${SEISSOL_SYSTEM_INCLUDES})

    target_include_directories(SeisSol-omp-offloaded PRIVATE ${SEISSOL_INCLUDES})
    target_include_directories(SeisSol-omp-offloaded SYSTEM PRIVATE ${SEISSOL_SYSTEM_INCLUDES})
    target_link_libraries(SeisSol-omp-offloaded PRIVATE easi::easi)

    install(TARGETS SeisSol-omp-offloaded DESTINATION ${CMAKE_INSTALL_PREFIX})
  endif()
  install(TARGETS SeisSol-device-lib DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(TARGETS device DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()


# Fortran compliler settings
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(SeisSol-lib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-cpp -ffree-line-length-none -fdefault-real-8 -Wno-unused-parameter>)
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # todo intel, is needed:  -align -align array64byte
  # todo -r8 -WB is needed for intel (8 byte precision for reals)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -r8 -WB")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    # NOTE:
    # 1. -Mpreprocess enables a C preprocessor macro
    # 2. -Mfree specifies fortran free format
    # 3. -r8 Interpret REAL variables as DOUBLE PRECISION
    # 4. BG (see, Initializer/preProcessorMacros.fpp) around logDebug, logInfo, etc.
  target_compile_options(SeisSol-lib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess -Mfree -r8 -DBG>)
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
  target_compile_options(SeisSol-lib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-cpp -ffree-form -fdefault-real-8 -DBG>)
endif()

# C++ compiler settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(SeisSol-lib PUBLIC -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas>)

  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  # todo remove std?, is ffreestanding needed?
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

  # Activate interprocedual optimization.
  #set_property(TARGET SeisSol-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE) 
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -pedantic -Wno-unused-parameter")

elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
    set(WARNINGS "--display_error_number --diag_suppress186")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "22.3.1")
      set(WARNINGS "${WARNINGS} --diag_suppress1")
    endif()

    # NOTE:
    # 1. --pending_instantiations=0 allows an infinite recursive template instantiation
    # 2. EIGEN_DONT_VECTORIZE=1 waiting for eigen3 support for nvhpc compiler collection w.r.t. vectorization
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wc,--pending_instantiations=0 ${WARNINGS} -DEIGEN_DONT_VECTORIZE=0")
endif()

find_package(FILESYSTEM REQUIRED)
target_link_libraries(SeisSol-lib PUBLIC std::filesystem)

# Generated code does only work without red-zone.
if (HAS_REDZONE)
  set_source_files_properties(
      ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone
  )
endif()

# adjust prefix name of executables
if ("${DEVICE_ARCH_STR}" STREQUAL "none")
  set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${HOST_ARCH_STR}_${ORDER}_${EQUATIONS}")
else()
  set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${DEVICE_ARCH_STR}_${DEVICE_BACKEND}_${ORDER}_${EQUATIONS}")
endif()


add_executable(SeisSol-bin src/main.cpp)
set_target_properties(SeisSol-bin PROPERTIES OUTPUT_NAME "SeisSol_${EXE_NAME_PREFIX}")
target_link_libraries(SeisSol-bin PUBLIC SeisSol-lib)

# SeisSol proxy-core
add_library(SeisSol-proxy-core auto_tuning/proxy/src/proxy_seissol.cpp)
target_link_libraries(SeisSol-proxy-core PUBLIC SeisSol-lib)

if (PROXY_PYBINDING)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    find_package(pybind11 2.6.2 REQUIRED)
  else()
    find_package(pybind11 2.3.0 REQUIRED)
  endif()

  pybind11_add_module(seissol_proxy_bindings auto_tuning/proxy/src/proxy_bindings.cpp)
  target_link_libraries(seissol_proxy_bindings PUBLIC SeisSol-proxy-core)
  file(COPY auto_tuning/proxy/src/proxy-runners DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()


if (WITH_GPU)
  if (CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
    target_link_options(SeisSol-bin PRIVATE "-mp=gpu")
  endif()
  if (OpenMP_Offloading_FOUND)
    target_link_libraries(SeisSol-bin PUBLIC SeisSol-omp-offloaded)
    target_link_libraries(SeisSol-proxy-core PUBLIC SeisSol-omp-offloaded)
  endif()
endif()

# C/C++ proxy interface
add_executable(SeisSol-proxy auto_tuning/proxy/src/proxy_main.cpp)
target_link_libraries(SeisSol-proxy PUBLIC SeisSol-proxy-core SeisSol-lib)
set_target_properties(SeisSol-proxy PROPERTIES OUTPUT_NAME "SeisSol_proxy_${EXE_NAME_PREFIX}")
install(TARGETS SeisSol-proxy DESTINATION ${CMAKE_INSTALL_PREFIX})

if (LIKWID)
  find_package(likwid REQUIRED)
  target_compile_definitions(SeisSol-proxy-core PUBLIC LIKWID_PERFMON)
  target_link_libraries(SeisSol-proxy-core PRIVATE likwid::likwid)
endif()

if (TESTING)
  enable_testing()
  include(cmake/doctest.cmake)

  # Coverage
  if(COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
    include(cmake/CodeCoverage.cmake)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(
            NAME SeisSol-coverage
            EXECUTABLE SeisSol-serial-test
            EXCLUDE "/usr/*"
                    "submodules/*"
                    "*/tests/*"
                    "external/*"
                    "*/yaml-cpp-install/*"
    )
  endif()

  set(seissol_test_sources
          src/tests/Model/TestModel.cpp
          src/tests/Initializer/TestInitializer.cpp
          src/tests/Numerical_aux/TestNumerical_aux.cpp
          src/tests/Geometry/TestGeometry.cpp
          src/tests/Kernel/TestKernel.cpp
          src/tests/SourceTerm/TestSourceTerm.cpp
          src/tests/Pipeline/TestPipeline.cpp
          src/tests/ResultWriter/TestResultWriter.cpp
          src/tests/DynamicRupture/TestDynamicRupture.cpp
          )


  if (TESTING_GENERATED)
    set(seissol_test_sources
            ${seissol_test_sources}
            ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/test-kernel.cpp
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

  if (NETCDF)
    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/source_loh.nrf
            ${CMAKE_CURRENT_BINARY_DIR}/Testing/source_loh.nrf
            COPYONLY
    )
    set(seissol_test_sources ${seissol_test_sources} ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/ReaderTest.cpp)
  endif()

  add_executable(SeisSol-serial-test
            ${seissol_test_sources}
            src/tests/test_main.cpp)
  target_link_libraries(SeisSol-serial-test PRIVATE SeisSol-lib)
  target_include_directories(SeisSol-serial-test PUBLIC external/)
  doctest_discover_tests(SeisSol-serial-test)

  # Avoid duplicate definition of FLOP counters
  target_compile_definitions(SeisSol-serial-test PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)
endif()

install(TARGETS SeisSol-bin DESTINATION ${CMAKE_INSTALL_PREFIX})
# https://blog.kitware.com/static-checks-with-cmake-cdash-iwyu-clang-tidy-lwyu-cpplint-and-cppcheck/
# https://ortogonal.github.io/cmake-clang-tidy/
