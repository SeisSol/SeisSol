cmake_minimum_required(VERSION 3.5)
project(SeisSol LANGUAGES C CXX Fortran)

message(WARNING "We do not provide support for CMake currently. Please use scons instead!")

# Library settings:
option(HDF5 "Use HDF5 library for data output" ON)
option(NETCDF "Use netcdf library for mesh input" ON)
option(METIS "Use metis for partitioning" ON)
option(MPI "Use MPI parallelization" ON)
option(OPENMP "Use OpenMP parallelization" ON)
option(ASAGI "Use asagi for material input" OFF)

# todo:
option(SIONLIB "Use sionlib for checkpointing" OFF)
option(MEMKIND "Use memkind library for hbw memory support" OFF)

# Other settings
set(ORDER "6" CACHE STRING "Convergence order.")
set_property(CACHE ORDER PROPERTY STRINGS "2" "3" "4" "5" "6" "7" "8")
set(NUMBER_OF_MECHANISMS "0" CACHE STRING "Number of mechanisms")
set(EQUATIONS "elastic" CACHE STRING "Equation set used")
set(ARCH "dsnb" CACHE STRING "Architecture")
set(DYNAMIC_RUPTURE_METHOD "quadrature" CACHE STRING "Dynamic rupture method.")
set_property(CACHE DYNAMIC_RUPTURE_METHOD PROPERTY STRINGS "quadrature" "cellaverage")
option(PLASTICITY "Use plasticity")
set(PLASTICITY_METHOD "nb" CACHE STRING "Dynamic rupture method.")
set_property(CACHE PLASTICITY_METHOD PROPERTY STRINGS nb ip)
set(NUMBER_OF_MULTIPLE_SIMULATIONS "1" CACHE STRING "Fuse multiple simulations in one run.")
option(COMMTHREAD "Use a communication thread for MPI+MP." OFF)
set(LOG_LEVEL "INFO" CACHE STRING "Log level for the code")
set_property(CACHE LOG_LEVEL PROPERTY STRINGS "DEBUG" "INFO" "WARNING" "ERROR")


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

# PDE-Settings
if ("${EQUATIONS}" STREQUAL "elastic" AND ${NUMBER_OF_MECHANISMS} GREATER 0)
  message(FATAL_ERROR "${EQUATIONS} does not support a NUMBER_OF_MECHANISMS > 0.")
endif()

MATH(EXPR NUMBER_OF_QUANTITIES "9 + 6 * ${NUMBER_OF_MECHANISMS}" )

include(src/sources.cmake)
add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/generated_code")

add_custom_command(
  COMMAND
  "/usr/bin/env" python3 ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/generate.py
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code
     "--arch" ${ARCH}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" "auto"
     "--multipleSimulations" ${NUMBER_OF_MULTIPLE_SIMULATIONS}
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

function(target_compiler_flags target)
  cmake_parse_arguments(
    TARGET_COMPILER_FLAGS
    ""
    "COMPILER"
    "PUBLIC;PRIVATE"
    ${ARGN}
    )
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL ${TARGET_COMPILER_FLAGS_COMPILER})
    target_compile_options(${target}
      PRIVATE ${TARGET_COMPILER_FLAGS_PRIVATE}
      PUBLIC ${TARGET_COMPILER_FLAGS_PUBLIC})
  endif()
endfunction()

# Set float size
if ("${ARCH}" MATCHES "^s.*")
  # single precision
  target_compile_definitions(SeisSol-lib PUBLIC REAL_SIZE=4)
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNED_REAL_SIZE=4)
elseif ("${ARCH}" MATCHES "^d.*")
  # double precision
  target_compile_definitions(SeisSol-lib PUBLIC REAL_SIZE=8)
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNED_REAL_SIZE=8)
else()
  message(FATAL_ERROR "Arch has to begin with s or d.")
endif()
  
# Arch
set(ALIGNMENT 0)
if ("${ARCH}" MATCHES "^(s|d)noarch$")
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNMENT=16)
elseif ("${ARCH}" MATCHES "^(s|d)wsm$")
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNMENT=16)
  target_compile_options(SeisSol-lib PUBLIC "-msse3")
elseif ("${ARCH}" MATCHES "^(s|d)snb$")
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNMENT=32)
  target_compile_options(SeisSol-lib PUBLIC "-mavx")
elseif ("${ARCH}" MATCHES "^(s|d)hsw$")
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNMENT=32)
  target_compiler_flags(SeisSol-lib COMPILER "Intel" PUBLIC "-xCORE-AVX2" "-fma")
  target_compiler_flags(SeisSol-lib COMPILER "GNU" PUBLIC "-mavx2" "-mfma")
elseif ("${ARCH}" MATCHES "^(s|d)knc$")
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNMENT=64)
  target_compile_options(SeisSol-lib PUBLIC "-mmic" "-fma")
elseif ("${ARCH}" MATCHES "^(s|d)knl$")
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNMENT=64)
  target_compiler_flags(SeisSol-lib COMPILER "Intel" PUBLIC "-xMIC-AVX512" "-fma")
  target_compiler_flags(SeisSol-lib COMPILER "GNU" PUBLIC "-mavx512f" "-mavx512cd" "-mavx512pf" "-mavx512er" "-mfma")
elseif ("${ARCH}" MATCHES "^(s|d)skx$")
  target_compile_definitions(SeisSol-lib PUBLIC ALIGNMENT=64)
  target_compiler_flags(SeisSol-lib COMPILER "Intel" PUBLIC "-xCORE-AVX512" "-fma")
  target_compiler_flags(SeisSol-lib COMPILER "GNU" PUBLIC "-march=skylake-avx512")
else ()
  message(FATAL_ERROR "Arch is unknown.")
  endif()

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
target_link_libraries(SeisSol-lib PUBLIC yaml-cpp
			      ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ImpalaJIT/libimpalajit.a 
			      )

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
  target_link_directories(SeisSol-lib PUBLIC ${PARMETIS_LIBRARY_DIR})
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

set_target_properties(SeisSol-bin PROPERTIES OUTPUT_NAME "SeisSol_${CMAKE_BUILD_TYPE}_${ARCH}_${EQUATIONS}_${ORDER}")

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
