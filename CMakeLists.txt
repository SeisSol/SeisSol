cmake_minimum_required(VERSION 3.9)

if(POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

project(SeisSol LANGUAGES C CXX Fortran)

# set hardware specific definition needed for seissol compilation
# 'process_users_input' returns the following:
#
#       switches: HDF5, NETCDF, METIS, MPI, OPENMP, ASAGI, SIONLIB, MEMKIND
#
#       user's input: HOST_ARCH, DEVICE_ARCH, DEVICE_SUB_ARCH,
#                     ORDER, NUMBER_OF_MECHANISMS, EQUATIONS,
#                     PRECISION, DYNAMIC_RUPTURE_METHOD,
#                     PLASTICITY, NUMBER_OF_FUSED_SIMULATIONS,
#                     MEMORY_LAYOUT, COMMTHREAD,
#                     LOG_LEVEL, LOG_LEVEL_MASTER,
#                     GEMM_TOOLS_LIST
#
#       derived input: REAL_SIZE_IN_BYTES, ALIGNMENT, ARCH_STRING, NUMBER_OF_QUANTITIES, WITH_GPU
#
include(cmake/process_users_input.cmake)
set(HARDWARE_DEFINITIONS "ALIGNMENT=${ALIGNMENT}"
                         "REAL_SIZE=${REAL_SIZE_IN_BYTES}"
                         "ALIGNED_REAL_SIZE=${REAL_SIZE_IN_BYTES}")

# TODO: Move one dir up.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

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

find_package(GemmTools REQUIRED)
include(src/sources.cmake)

if(GemmTools_INCLUDE_DIRS)
  target_include_directories(SeisSol-lib PUBLIC ${GemmTools_INCLUDE_DIRS})
endif()
if(GemmTools_LIBRARIES)
  target_link_libraries(SeisSol-lib PUBLIC ${GemmTools_LIBRARIES})
endif()
if(GemmTools_COMPILER_DEFINITIONS)  
  target_compile_definitions(SeisSol-lib PUBLIC ${GemmTools_COMPILER_DEFINITIONS})
endif()

add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/generated_code")

add_custom_command(
  COMMAND
  "/usr/bin/env" python3 ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/generate.py
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code
     "--host_arch" ${HOST_ARCH_STR}
     "--device_arch" ${DEVICE_ARCH_STR}
     "--device_sub_arch" ${DEVICE_SUB_ARCH}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" ${MEMORY_LAYOUT}
     "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
     "--dynamicRuptureMethod" ${DYNAMIC_RUPTURE_METHOD}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
     "--gemm_tools" ${GEMM_TOOLS_LIST}
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
       generated_code/anisotropic.py
       generated_code/SurfaceDisplacement.py
       generated_code/NodalBoundaryConditions.py
    OUTPUT src/generated_code/subroutine.h
       src/generated_code/tensor.cpp
       src/generated_code/subroutine.cpp
       src/generated_code/gpulike_subroutine.cpp
       src/generated_code/tensor.h
       src/generated_code/init.cpp
       src/generated_code/init.h
       src/generated_code/kernel.h
       src/generated_code/kernel.cpp
       COMMENT "Codegen for tensor stuff."
       )

if (PLASTICITY)
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY)
  set(PLASTICITY_NAME_SUFFIX "_plasticity")
endif()

if (PLASTICITY_METHOD STREQUAL "ip")
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY_IP)
elseif (PLASTICITY_METHOD STREQUAL "nb")
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY_NB)
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

# enable interproc. opts for small cores
#if cpu in ['knc', 'knl', 'skx']:
#    flags.extend(['-ip'])

# Libs
include(ExternalProject)

find_package(YAML-CPP 0.5.3 QUIET)
if (YAML-CPP_FOUND)
  target_link_libraries(SeisSol-lib PUBLIC ${YAML_CPP_LIBRARIES})
  target_include_directories(SeisSol-lib PUBLIC ${YAML_CPP_INCLUDE_DIR})
else()
  ExternalProject_Add(yaml-cpp-build
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/submodules/yaml-cpp/
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-install
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DYAML_CPP_BUILD_TOOLS=OFF
    -DYAML_CPP_BUILD_TESTS=OFF
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  )
  add_dependencies(SeisSol-lib yaml-cpp-build)
  target_link_libraries(SeisSol-lib PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-install/lib/libyaml-cpp.a")
  target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-install/include/)
endif()


ExternalProject_Add(ImpalaJIT-build
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/submodules/ImpalaJIT/
  INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT-install
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)
add_dependencies(SeisSol-lib ImpalaJIT-build)
target_link_libraries(SeisSol-lib PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT-install/lib/libimpalajit.a")
target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT-install/include/)

# todo: if parallelization == omp/hybrid
if (OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(SeisSol-lib PUBLIC OpenMP::OpenMP_CXX)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  target_compile_definitions(SeisSol-lib PUBLIC OMP OMPI_SKIP_MPICXX)
endif()

if (MPI)
  find_package(MPI REQUIRED)

  target_include_directories(SeisSol-lib SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC MPI::MPI_C)

  target_include_directories(SeisSol-lib SYSTEM PUBLIC ${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC MPI::MPI_Fortran)

  target_compile_definitions(SeisSol-lib PUBLIC USE_MPI PARALLEL)
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
  pkg_check_modules(ASAGI REQUIRED asagi) # asagi_nompi?
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

if (WITH_GPU)
  add_subdirectory(submodules/Device)
  target_include_directories(SeisSol-lib PUBLIC submodules/Device)
  target_link_libraries(SeisSol-lib PUBLIC device)
  target_compile_definitions(SeisSol-lib PUBLIC ACL_DEVICE)
endif()

target_include_directories(SeisSol-lib PUBLIC
				   ${CMAKE_CURRENT_SOURCE_DIR}/src
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/async
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen3
				   # todo: use reexport from easi?
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/easi/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
				   ${CMAKE_CURRENT_BINARY_DIR}/src/
)

target_include_directories(SeisSol-lib SYSTEM PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen3
)

target_compile_definitions(SeisSol-lib PUBLIC
				   CONVERGENCE_ORDER=${ORDER}
				   NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
				   NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
				   ENABLE_MATRIX_PREFETCH
				   )

# Fortran compliler settings
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(SeisSol-lib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none -fdefault-real-8 -fopenmp -Wno-unused-parameter>)
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # todo intel, is needed:  -align -align array64byte
  # todo -r8 -WB is needed for intel (8 byte precision for reals)
  # todo openmp only if desired
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -r8 -WB -qopenmp")
endif()

# C++ compiler settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(SeisSol-lib PUBLIC -fopenmp -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas>)
  target_link_libraries(SeisSol-lib PUBLIC "-fopenmp")

  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  # todo remove std?, is ffreestanding needed?
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -qopenmp")

  # only if openmp
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -qopenmp")

  # Activate interprocedual optimization.
  #set_property(TARGET SeisSol-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE) 
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp=libomp -Wall -Wextra -pedantic")
endif()

# Generated code does only work without red-zone.
if (HAS_REDZONE)
  set_source_files_properties(src/generated_code/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone)
endif()

# adjust prefix name of executables
if ("${DEVICE_ARCH_STR}" STREQUAL "none")
  set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${HOST_ARCH_STR}_${ORDER}_${EQUATIONS}")
else()
  set(EXE_NAME_PREFIX "${CMAKE_BUILD_TYPE}_${DEVICE_ARCH_STR}_${ORDER}_${EQUATIONS}")
endif()


add_executable(SeisSol-bin src/main.cpp)
target_link_libraries(SeisSol-bin PUBLIC SeisSol-lib)
set_target_properties(SeisSol-bin PROPERTIES OUTPUT_NAME "SeisSol_${EXE_NAME_PREFIX}${PLASTICITY_NAME_SUFFIX}")

add_executable(SeisSol-proxy
        auto_tuning/proxy/src/proxy_seissol.cpp
        auto_tuning/proxy/src/flop_counter.cpp
)
target_link_libraries(SeisSol-proxy PUBLIC SeisSol-lib)
set_target_properties(SeisSol-proxy PROPERTIES OUTPUT_NAME "SeisSol_proxy_${EXE_NAME_PREFIX}")


if (TESTING)
  include(cmake/testing.cmake)

  # Testing
  set(CXXTEST_TESTGEN_ARGS "--template=${CMAKE_CURRENT_SOURCE_DIR}/src/tests/mpirunner.tpl")
  find_package(CxxTest REQUIRED)

  # maybe enable CXXTEST_HAVE_EH CXXTEST_ABORT_ON_FAIL
  enable_testing()

  # what about the unit tests for the lts layout?

  if (TESTING_GENERATED)
    CXXTEST_ADD_TEST_MPI(
            test_kernel_test_suite
            1
            test_kernel.cpp
            ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/KernelTest.t.h
    )

    target_link_libraries(test_kernel_test_suite PRIVATE SeisSol-lib)
    target_include_directories(test_kernel_test_suite PRIVATE ${CXXTEST_INCLUDE_DIR})
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

  if (NETCDF)
  configure_file(
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/source_loh.nrf 
          ${CMAKE_CURRENT_BINARY_DIR}/Testing/source_loh.nrf 
          COPYONLY
  )
  set(SeisSol_NETCDF_TEST_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Reader/NRFReader.t.h)
  endif()

  CXXTEST_ADD_TEST_MPI(
          test_serial_test_suite
          1
          test_serial.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Numerical_aux/Functions.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Numerical_aux/Quadrature.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Numerical_aux/Transformations.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Physics/PointSource.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Model/GodunovState.t.h
	      ${SeisSol_NETCDF_TEST_FILES}
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Geometry/MeshRefiner.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Geometry/VariableSubsampler.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Geometry/TriangleRefiner.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Initializer/time_stepping/LTSWeights.t.h
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Initializer/PointMapper.t.h
  )
  target_link_libraries(test_serial_test_suite PRIVATE SeisSol-lib)
  target_include_directories(test_serial_test_suite PRIVATE ${CXXTEST_INCLUDE_DIR})

  CXXTEST_ADD_TEST_MPI(
          test_parallel_test_suite
          4
          test_parallel.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/minimal/Minimal.t.h
  )
  target_link_libraries(test_parallel_test_suite PRIVATE SeisSol-lib)
  target_include_directories(test_parallel_test_suite PRIVATE ${CXXTEST_INCLUDE_DIR})
  target_compile_definitions(test_parallel_test_suite PRIVATE
    SEISSOL_TESTS=\"${CMAKE_CURRENT_SOURCE_DIR}/src/tests/\"
    )
endif()

# https://blog.kitware.com/static-checks-with-cmake-cdash-iwyu-clang-tidy-lwyu-cpplint-and-cppcheck/
# https://ortogonal.github.io/cmake-clang-tidy/
