cmake_minimum_required(VERSION 3.5)
project(SeisSol LANGUAGES C CXX Fortran)

option(HDF5 "Use HDF5 library for data output" ON)
option(NETCDF "Use netcdf library for mesh input" ON)
option(METIS "Use metis for partitioning" ON)
option(MPI "Use MPI parallelization" ON)
option(OPENMP "Use OpenMP parallelization" ON)
# todo:
option(ASAGI "Use asagi for material input" OFF) # needed
option(SIONLIB "Use sionlib for checkpointing" OFF)
option(MEMKIND "Use memkind library for hbw memory support" OFF)
option(COMMTHREAD "Use a communication thread for MPI+MP." OFF)

set(ORDER 6 CACHE INT "Convergence order.")
set_property(CACHE ORDER PROPERTY STRINGS 2 3 4 5 6 7 8)
set(NUMBER_OF_MECHANISMS 0 CACHE INT "Number of mechanisms")
set(EQUATIONS "elastic" CACHE STRING "Equation set used")
set(ARCH "dsnb" CACHE STRING "Architecture")
set(DYNAMIC_RUPTURE_METHOD "quadrature" CACHE STRING "Dynamic rupture method.")
set_property(CACHE DYNAMIC_RUPTURE_METHOD PROPERTY STRINGS "quadrature" "cellaverage")
set(PLASTICITY_METHOD "nb" CACHE STRING "Dynamic rupture method.")
set_property(CACHE PLASTICITY_METHOD PROPERTY STRINGS nb ip)
set(NUMBER_OF_MULTIPLE_SIMULATIONS 1 CACHE INT "Fuse multiple simulations in one run.")

# TODO: Set number of quantities depending on equations

# TODO: Move one dir up.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
#set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/cmake_build)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
if ("${EQUATIONS}" MATCHES "elastic" AND ${NUMBER_OF_MECHANISMS} GREATER 0)
  message(FATAL_ERROR "${EQUATIONS} does not support a NUMBER_OF_MECHANISMS > 0.")
endif()

MATH(EXPR NUMBER_OF_QUANTITIES "9 + 6 * ${NUMBER_OF_MECHANISMS}" )

# Set float size
if ("${ARCH}" MATCHES "^s.*")
  # single precision
  set(REAL_SIZE 4)
  set(ALIGNED_REAL_SIZE 4)
elseif ("${ARCH}" MATCHES "^d.*")
  # double precision
  set(REAL_SIZE 8)
  set(ALIGNED_REAL_SIZE 8)
else()
  message(FATAL_ERROR "Arch has to begin with s or d.")
endif()
  
# Arch
set(ALIGNMENT 0)
if ("${ARCH}" MATCHES "^(s|d)noarch$")
  set(ALIGMENT 16)
elseif ("${ARCH}" MATCHES "^(s|d)wsm$")
  set(ALIGNMENT 16)
  # -msse3
elseif ("${ARCH}" MATCHES "^(s|d)snb$")
  set(ALIGNMENT 32)
  # -mavx
elseif ("${ARCH}" MATCHES "^(s|d)hsw$")
  set(ALIGNMENT 32)
  # intel:
  # intel -> -xCore-AVX2 -fma
  # gcc -> -mavx2 -mfma
elseif ("${ARCH}" MATCHES "^(s|d)knc$")
  set(ALIGNMENT 64)
  # -mmic -fma
elseif ("${ARCH}" MATCHES "^(s|d)knl$")
  set(ALIGNMENT 64)
  # intel -> -xMIC-AVX512 -fma
  # gcc -> -mavx512f -mavx512cd -mavx512pf -mavx512er -mfma
elseif ("${ARCH}" MATCHES "^(s|d)skx$")
  set(ALIGNMENT 64)
  # intel -> -xCORE-AVX512 -fma
  # gcc -> -march=skylake-avx512
else ()
  message(FATAL_ERROR "Arch is unknown.")
  endif()

# enable interproc. opts for small cores
#if cpu in ['knc', 'knl', 'skx']:
#    flags.extend(['-ip'])


# todo: guess memlayout
# todo create subdir $BIN/src/generated_code
# todo compile subroutine with -mno-red-zone
add_custom_command(
  COMMAND
  "/usr/bin/env" python3 ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/generate.py
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code
     "--arch" ${ARCH}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" /work_fast/krenz/src/SeisSol-CMake/auto_tuning/config/dense.xml
     "--multipleSimulations" ${NUMBER_OF_MULTIPLE_SIMULATIONS}
     "--dynamicRuptureMethod" ${DYNAMIC_RUPTURE_METHOD}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/generated_code
     DEPENDS
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
    OUTPUT
       src/generated_code/subroutine.h
       src/generated_code/tensor.cpp
       src/generated_code/subroutine.cpp
       src/generated_code/tensor.h
       src/generated_code/init.cpp
       src/generated_code/init.h
       src/generated_code/kernel.h
       src/generated_code/kernel.cpp
       COMMENT "Codegen for tensor stuff."
       )


add_library(SeisSol-lib
src/Initializer/ParameterDB.cpp
src/Initializer/PointMapper.cpp
src/Initializer/GlobalData.cpp
src/Initializer/InternalState.cpp
src/Initializer/MemoryAllocator.cpp
src/Initializer/CellLocalMatrices.cpp

# TODO: Only if MPI?
src/Initializer/time_stepping/LtsWeights.cpp

src/Initializer/time_stepping/LtsLayout.cpp
src/Initializer/tree/Lut.cpp
src/Initializer/MemoryManager.cpp
src/Initializer/InitialFieldProjection.cpp
src/Modules/Modules.cpp
src/Modules/ModulesC.cpp
src/Model/common.cpp
src/Numerical_aux/Transformation.cpp

src/generated_code/subroutine.h
src/generated_code/tensor.cpp
src/generated_code/subroutine.cpp
src/generated_code/tensor.h
src/generated_code/init.cpp
#src/generated_code/KernelTest.t.h
src/generated_code/init.h
src/generated_code/kernel.h
src/generated_code/kernel.cpp
src/
# Equations/viscoelastic2/Model/Setup.cpp
# Equations/viscoelastic2/Kernels/Neighbor.cpp
# Equations/viscoelastic2/Kernels/Local.cpp
# Equations/viscoelastic2/Kernels/Time.cpp
src/Equations/elastic/Model/Setup.cpp
src/Equations/elastic/Kernels/Neighbor.cpp
src/Equations/elastic/Kernels/Local.cpp
src/Equations/elastic/Kernels/Time.cpp
src/Equations/elastic/Physics/InitialField.cpp
# Equations/viscoelastic/Model/Setup.cpp
# Equations/viscoelastic/Kernels/Neighbor.cpp
# Equations/viscoelastic/Kernels/Local.cpp
# Equations/viscoelastic/Kernels/Time.cpp
src/Solver/Simulator.cpp
src/Solver/FreeSurfaceIntegrator.cpp
src/Solver/Interoperability.cpp
src/Solver/time_stepping/MiniSeisSol.cpp
src/Solver/time_stepping/TimeCluster.cpp
src/Solver/time_stepping/TimeManager.cpp
src/Kernels/DynamicRupture.cpp
src/Kernels/Plasticity.cpp
src/Kernels/TimeCommon.cpp
src/Kernels/Receiver.cpp
src/SeisSol.cpp
src/SourceTerm/Manager.cpp

src/SourceTerm/NRFReader.cpp # if netCDF

src/SourceTerm/PointSource.cpp
src/Parallel/Pin.cpp
src/Parallel/MPI.cpp
src/Parallel/mpiC.cpp
src/Parallel/FaultMPI.cpp
src/Geometry/GambitReader.cpp

# TODO: Only if hdf5?
src/Geometry/PUMLReader.cpp

src/Geometry/MeshReaderFBinding.cpp
src/Geometry/MeshTools.cpp
src/Monitoring/FlopCounter.cpp
src/Monitoring/LoopStatistics.cpp
src/Reader/readparC.cpp
#Reader/StressReaderC.cpp
#Reader/AsagiModule.cpp
src/Checkpoint/Manager.cpp

# TODO: Only if mpi?
src/Checkpoint/mpio/Wavefield.cpp
src/Checkpoint/mpio/FaultAsync.cpp
src/Checkpoint/mpio/Fault.cpp
src/Checkpoint/mpio/WavefieldAsync.cpp

# Checkpoint/sionlib/Wavefield.cpp
# Checkpoint/sionlib/Fault.cpp

# TODO: Only if HDF5!
src/Checkpoint/h5/Wavefield.cpp
src/Checkpoint/h5/Fault.cpp

src/Checkpoint/Backend.cpp
src/Checkpoint/Fault.cpp
src/Checkpoint/posix/Wavefield.cpp
src/Checkpoint/posix/Fault.cpp
src/ResultWriter/AnalysisWriter.cpp
src/ResultWriter/FreeSurfaceWriterExecutor.cpp
src/ResultWriter/PostProcessor.cpp
src/ResultWriter/FaultWriterC.cpp
src/ResultWriter/ReceiverWriter.cpp
src/ResultWriter/FaultWriterExecutor.cpp
src/ResultWriter/FaultWriter.cpp
src/ResultWriter/WaveFieldWriter.cpp
src/ResultWriter/FreeSurfaceWriter.cpp

# Fortran:
src/Monitoring/bindMonitoring.f90
src/Geometry/mpiextractmesh.f90
src/Geometry/allocate_mesh.f90
src/Geometry/MeshReaderCBinding.f90
src/Solver/close_seissol.f90
src/Solver/calc_deltat.f90
src/Solver/mpiexchangevalues.f90
src/Solver/prak_clif_mod.f90
src/Solver/calc_seissol.f90
src/Solver/f_ctof_bind_interoperability.f90
src/Solver/plasticity.f90
src/Solver/f_ftoc_bind_interoperability.f90
src/Numerical_aux/quadpoints.f90
src/Numerical_aux/jacobinormal.f90
src/Numerical_aux/convertxieta2xy.f90
src/Numerical_aux/create_fault_rotationmatrix.f90
src/Numerical_aux/trilinearinterpolation.f90
src/Numerical_aux/typesdef.f90
src/Numerical_aux/dgbasis.f90
src/Numerical_aux/gauss.f90
src/Numerical_aux/operators.f90
src/Modules/ModulesF.f90
src/seissolxx.f90
src/Physics/ini_model.f90
src/Physics/Evaluate_friction_law.f90
src/Physics/ini_model_DR.f90
src/Reader/readpar.f90
src/Reader/read_backgroundstress.f90
src/ResultWriter/inioutput_seissol.f90
src/ResultWriter/magnitude_output.f90
src/ResultWriter/output_rupturefront.f90
src/ResultWriter/ini_faultoutput.f90
src/ResultWriter/energies.f90
src/ResultWriter/FaultWriterF.f90
src/ResultWriter/faultoutput.f90
src/ResultWriter/common_fault_receiver.f90
src/ResultWriter/receiver.f90
src/Initializer/dg_setup.f90
src/Initializer/ini_optionalfields.f90
src/Initializer/ini_seissol.f90
src/Parallel/mpiF.f90
)

#add_dependencies(SeisSol-lib yateto-codegen)


# Libs
find_package(yaml-cpp REQUIRED)

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
  include(FindMPI)
  target_include_directories(SeisSol-lib PUBLIC ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC ${MPI_CXX_LIBRARIES})
  target_link_libraries(SeisSol-lib PUBLIC ${MPI_C_LIBRARIES})
  target_include_directories(SeisSol-lib PUBLIC ${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(SeisSol-lib PUBLIC ${MPI_Fortran_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_MPI PARALLEL)
endif()

# TODO: HDF5 should be optional?
# TODO: Check constraints, seee
# https://github.com/Kitware/CMake/blob/master/Modules/FindHDF5.cmake
if (HDF5)
  find_package(HDF5 REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC ${HDF5_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_HDF)
endif()  

if (NETCDF)
  find_package(NetCDF REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC NetCDF::NetCDF)
  target_compile_definitions(SeisSol-lib PUBLIC USE_NETCDF)
endif()

# Parmetis
if (METIS)
  find_package(ParMETIS REQUIRED)
  target_include_directories(SeisSol-lib PUBLIC ${PARMETIS_INCLUDE_DIRS})
  target_link_libraries(SeisSol-lib PUBLIC ${PARMETIS_LIBRARIES})
  target_compile_definitions(SeisSol-lib PUBLIC USE_METIS)
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
				   # todo: Switch kernels?
				   ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic
				 
)

target_compile_definitions(SeisSol-lib PUBLIC
				   REAL_SIZE=${REAL_SIZE}
				   ALIGNMENT=${ALIGNMENT}
				   CONVERGENCE_ORDER=${ORDER}
				   NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
				   NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
				   # TODO: Remove next one
				   VERSION_STRING=""
				   __SSE3__
				   __AVX__
				   __AVX2__
				   GLM_FORCE_CXX98
				   ENABLE_MATRIX_PREFETCH
				   LOGLEVEL=3
				   LOGLEVEL0=3
				   # Are the following needed?
				   GENERATED_KERNELS
				   GENERATEDKERNELS
				   CLUSTERED_LTS
				   HOST_NAME_MAX=100 # wtf?

				   # some settings
				   USE_PLASTIC_NB
				   )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  # todo remove std?
  # todo is ffreestanding needed?
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  # todo intel, is needed:  -align -align array64byte
  # todo -r8 -WB is needed for intel (8 byte precision for reals)
  # todo openmp only if desired
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -r8 -WB -qopenmp")
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -qopenmp")

  # Activate interprocedual optimization.
  #set_property(TARGET SeisSol-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE) 
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "CRAY")
  # TODO: Cray
endif()


# Generated code does only work without red-zone.
set_source_files_properties(src/generated_code/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone)

add_executable(SeisSol-bin src/main.cpp)
target_link_libraries(SeisSol-bin PUBLIC SeisSol-lib)

# Testing
find_package(CxxTest REQUIRED)

if(CXXTEST_FOUND)
  # todo: if mpi set CXXTEST_COMMAND to mpirun -np 1 %t
  # maybe enable CXXTEST_HAVE_EH CXXTEST_ABORT_ON_FAIL
  enable_testing()
  set(CXXTEST_USE_PYTHON TRUE)

  # what about the unit tests for the lts layout?
  CXXTEST_ADD_TEST(test_generated_kernels
    test_generated_kernels.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/KernelTest.t.h)

  target_link_libraries(test_generated_kernels PRIVATE SeisSol-lib) 
  target_include_directories(test_generated_kernels PRIVATE ${CXXTEST_INCLUDE_DIR})

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
