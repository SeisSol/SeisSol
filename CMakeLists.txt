cmake_minimum_required(VERSION 3.5)
project(SeisSol LANGUAGES C CXX Fortran)
message(WARNING "We do not provide support for CMake currently. Please use scons instead!")


# set hardware specific definition needed for seissol for compilation
# 'process_users_input' returns the following:
#
#       switches: HDF5, NETCDF, METIS, MPI, OPENMP, ASAGI, SIONLIB, MEMKIND
#
#       user's input: ORDER, NUMBER_OF_MECHANISMS, EQUATIONS,
#                     ARCH, PRECISION, DYNAMIC_RUPTURE_METHOD,
#                     PLASTICITY, NUMBER_OF_FUSED_SIMULATIONS,
#                     MEMEORY_LAYOUT, COMMTHREAD, LOG_LEVEL,
#                     ACCELERATOR_TYPE
#
#       derived input: REAL_SIZE_IN_BYTES, ALIGNMENT, ARCH_STRING, NUMBER_OF_QUANTITIES
#
include(cmake/process_users_input.cmake)
set(HARDWARE_SPECIFIC_DEFINITIONS "ALIGNMENT=${ALIGNMENT}"
                                  "REAL_SIZE=${REAL_SIZE_IN_BYTES}"
                                  "ALIGNED_REAL_SIZE=${REAL_SIZE_IN_BYTES}")


# TODO: Move one dir up.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
  message(STATUS "Set build type to Debug a none was supplied.")
endif()


find_package(Git QUIET)
# Generate version.h
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE PACKAGE_GIT_VERSION
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file (
  "${PROJECT_SOURCE_DIR}/src/version.h"
  "${PROJECT_BINARY_DIR}/src/version.h"
  )

include(src/sources.cmake)
add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/generated_code")

add_custom_command(
  COMMAND
  "/usr/bin/env" python3 ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/generate.py
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code
     "--arch" ${ARCH_STRING}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" ${MEMEORY_LAYOUT}
     "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
     "--dynamicRuptureMethod" ${DYNAMIC_RUPTURE_METHOD}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
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
       generated_code/SurfaceDisplacement.py
    OUTPUT src/generated_code/subroutine.h
       src/generated_code/tensor.cpp
       src/generated_code/subroutine.cpp
       src/generated_code/tensor.h
       src/generated_code/init.cpp
       src/generated_code/init.h
       src/generated_code/kernel.h
       src/generated_code/kernel.cpp
       COMMENT "Codegen for tensor stuff."
       )


if (PLASTICITY)
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY)
endif()

if ("${PlasticityMethod}" STREQUAL "ip")
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY_IP)
elseif ("${PlasticityMethod}" STREQUAL "nb")
  target_compile_definitions(SeisSol-lib PUBLIC USE_PLASTICITY_NB)
endif()


# Find appropriate compiler flags based on the target computer architecture
# and the target compiler
# 'get_arch_specific_flags' returns: CPU_ARCH_SPECIFIC_FLAGS
include(cmake/cpu_arch_flags.cmake)
get_arch_specific_flags(${ARCH} ${CMAKE_CXX_COMPILER_ID})


# set hardware/compiler specific definitions and flags
target_compile_definitions(SeisSol-lib PUBLIC ${HARDWARE_SPECIFIC_DEFINITIONS})
target_compile_options(SeisSol-lib PUBLIC ${CPU_ARCH_SPECIFIC_FLAGS})


if (LOG_LEVEL STREQUAL "DEBUG")
  target_compile_definitions(SeisSol-lib PUBLIC LOGLEVEL=3 LOG_LEVEL=3 LOGLEVEL0=3)
elseif (LOG_LEVEL STREQUAL "INFO")
  target_compile_definitions(SeisSol-lib PUBLIC LOGLEVEL=2 LOG_LEVEL=2 LOGLEVEL0=2)
elseif (LOG_LEVEL STREQUAL "WARNING")
  target_compile_definitions(SeisSol-lib PUBLIC LOGLEVEL=1 LOG_LEVEL=1 LOGLEVEL0=1)
elseif (LOG_LEVEL STREQUAL "ERROR")
  target_compile_definitions(SeisSol-lib PUBLIC LOGLEVEL=0 LOG_LEVEL=0 LOGLEVEL0=0)
endif()
  
#set_property(CACHE LOG_LEVEL PRPERTY STRINGS "DEBUG" "INFO" "WARNING" "ERROR")
# enable interproc. opts for small cores
#if cpu in ['knc', 'knl', 'skx']:
#    flags.extend(['-ip'])

# Libs
find_package(yaml-cpp REQUIRED)

find_package(Eigen3) # maybe enfore 3.3?
if (EIGEN3_FOUND)
  target_include_directories(SeisSol-lib PUBLIC ${EIGEN3_INCLUDE_DIR})
  target_compile_definitions(SeisSol-lib PUBLIC HAS_EIGEN3)
else()
endif()

# TODO: Write export for impalaJIT
#target_link_libraries(SeisSol-lib PUBLIC yaml-cpp
				  #${CMAKE_CURRENT_SOURCE_DIR}/submodules/ImpalaJIT/libimpalajit.a 
				  #)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/ImpalaJIT)
target_link_libraries(SeisSol-lib PUBLIC yaml-cpp impalajit-static)



# todo: if parallelization == omp/hybrid
if (OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(SeisSol-lib PUBLIC OpenMP::OpenMP_CXX)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  target_compile_definitions(SeisSol-lib PUBLIC OMP OMPI_SKIP_MPICXX)
endif()

if (MPI)
  find_package(MPI REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC ${MPI_CXX_LIBRARIES})
  target_link_libraries(SeisSol-lib PUBLIC ${MPI_C_LIBRARIES})
  target_include_directories(SeisSol-lib PUBLIC ${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC ${MPI_Fortran_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_MPI PARALLEL)
endif()

if (COMMTHREAD)
  target_compile_definitions(SeisSol-lib PUBLIC USE_COMM_THREAD)
endif()

#set(HDF5_PREFER_PARALLEL True)
if (NETCDF)
  find_package(NetCDF REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${NetCDF_INCLUDE_DIRS})
  # todo: HL destroys search for hdf5 locally, set flags manually
  target_link_libraries(SeisSol-lib PUBLIC ${NetCDF_LIBRARY} "-lhdf5" "-lhdf5_hl")
  target_compile_definitions(SeisSol-lib PUBLIC USE_NETCDF)
endif()

if (HDF5)
  if (MPI)
    set(HDF5_PREFER_PARALLEL True)
  endif()
  find_package(HDF5 REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_HDF)
endif()  

# Parmetis
if (METIS)
  find_package(ParMETIS REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${PARMETIS_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC ${PARMETIS_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_METIS)
endif()

if (ASAGI)
  # todo warn if netcdf is off
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(ASAGI REQUIRED asagi) # asagi_nompi?
  target_compile_definitions(SeisSol-lib PUBLIC USE_ASAGI)
  target_link_libraries(SeisSol-lib PUBLIC ${ASAGI_STATIC_LDFLAGS})
  target_include_directories(SeisSol-lib PUBLIC ${ASAGI_INCLUDE_DIRS})
  target_compile_options(SeisSol-lib PUBLIC ${ASAGI_CFLAGS} ${ASAGI_CFLAGS_OTHER})
endif()

target_include_directories(SeisSol-lib PUBLIC
				   ${CMAKE_CURRENT_SOURCE_DIR}/src
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/async
				   # todo: use reexport from easi?
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/easi/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yaml-cpp/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ImpalaJIT/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/submodules/glm
				   ${CMAKE_CURRENT_BINARY_DIR}/src/
)

target_compile_definitions(SeisSol-lib PUBLIC
				   CONVERGENCE_ORDER=${ORDER}
				   NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
				   NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
				   GLM_FORCE_CXX98
				   ENABLE_MATRIX_PREFETCH
				   LOGLEVEL=3
				   LOGLEVEL0=3
				   # Are the following needed?
				   GENERATED_KERNELS
				   GENERATEDKERNELS
				   CLUSTERED_LTS
				   HOST_NAME_MAX=100 # wtf?
				   )

# Fortran compliler settings				 
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fdefault-real-8 -fopenmp")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # todo intel, is needed:  -align -align array64byte
  # todo -r8 -WB is needed for intel (8 byte precision for reals)
  # todo openmp only if desired
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -r8 -WB -qopenmp")
endif()

# C++ compiler settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(SeisSol-lib PUBLIC -fopenmp -Wall -Wextra -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wno-unused-parameter -Wno-unknown-pragmas>)
  target_link_options(SeisSol-lib PUBLIC -fopenmp)

  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  # todo remove std?, is ffreestanding needed?
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

  # only if openmp
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -qopenmp")

  # Activate interprocedual optimization.
  #set_property(TARGET SeisSol-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE) 
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp=libomp")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "CRAY")
  # TODO: Cray
endif()


# Generated code does only work without red-zone.
set_source_files_properties(src/generated_code/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone)

add_executable(SeisSol-bin src/main.cpp)
target_link_libraries(SeisSol-bin PUBLIC SeisSol-lib)

set_target_properties(SeisSol-bin PROPERTIES OUTPUT_NAME "SeisSol_${CMAKE_BUILD_TYPE}_${ARCH_STRING}_${EQUATIONS}_${ORDER}")

# Testing
find_package(CxxTest)

if(CXXTEST_FOUND)
  # todo: if mpi set CXXTEST_COMMAND to mpirun -np 1 %t
  # maybe enable CXXTEST_HAVE_EH CXXTEST_ABORT_ON_FAIL
  enable_testing()
  set(CXXTEST_USE_PYTHON TRUE)

  # what about the unit tests for the lts layout?

  # test is generated...
  #CXXTEST_ADD_TEST(test_generated_kernels
  #  test_generated_kernels.cpp
  #  ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/KernelTest.t.h)
  #target_link_libraries(test_generated_kernels PRIVATE SeisSol-lib) 
  #target_include_directories(test_generated_kernels PRIVATE ${CXXTEST_INCLUDE_DIR})

  CXXTEST_ADD_TEST(test_serial_test_suite
    test_serial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Numerical_aux/Functions.t.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Numerical_aux/Quadrature.t.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/Physics/PointSource.t.h
    )
  target_link_libraries(test_serial_test_suite PRIVATE SeisSol-lib)
  target_include_directories(test_serial_test_suite PRIVATE ${CXXTEST_INCLUDE_DIR})

  CXXTEST_ADD_TEST(test_parallel_test_suite
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
