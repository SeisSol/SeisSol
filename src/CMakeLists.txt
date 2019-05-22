cmake_minimum_required(VERSION 3.5)
project(SeisSol LANGUAGES C CXX Fortran)

option(hdf5 "Use HDF5 library for data output" ON)
option(netcdf "Use netcdf library for mesh input" ON)
option(metis "Use metis for partitioning" ON)
# todo:
option(asagi "Use asagi for material input" OFF) # needed
option(sionlib "Use sionlib for checkpointing" OFF)
option(memkind "Use memkind library for hbw memory support" OFF)

# TODO: Move one dir up.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake" ${CMAKE_MODULE_PATH})
#set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
#set (LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# todo remove std?
# todo is ffreestanding needed?
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ffreestanding")

    #enable_language(Fortran)
# todo intel, is needed:  -align -align array64byte
# todo -r8 -WB is needed for intel (8 byte precision for reals)
# todo openmp only if desired
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -r8 -WB -qopenmp")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -qopenmp")


find_package(Git QUIET)
# Generate version.h
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE PACKAGE_GIT_VERSION
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file (
  "${PROJECT_SOURCE_DIR}/version.h"
  "${PROJECT_BINARY_DIR}/version.h"
  )


# TODO: Generate code

add_executable(SeisSol
Initializer/ParameterDB.cpp
Initializer/PointMapper.cpp
Initializer/GlobalData.cpp
Initializer/InternalState.cpp
Initializer/MemoryAllocator.cpp
Initializer/CellLocalMatrices.cpp

# TODO: Only if MPI?
Initializer/time_stepping/LtsWeights.cpp

Initializer/time_stepping/LtsLayout.cpp
Initializer/tree/Lut.cpp
Initializer/MemoryManager.cpp
Initializer/InitialFieldProjection.cpp
Modules/Modules.cpp
Modules/ModulesC.cpp
Model/common.cpp
Numerical_aux/Transformation.cpp

generated_code/subroutine.h
generated_code/tensor.cpp
generated_code/subroutine.cpp
generated_code/tensor.h
generated_code/init.cpp
generated_code/KernelTest.t.h
generated_code/init.h
generated_code/kernel.h
generated_code/kernel.cpp

# Equations/viscoelastic2/Model/Setup.cpp
# Equations/viscoelastic2/Kernels/Neighbor.cpp
# Equations/viscoelastic2/Kernels/Local.cpp
# Equations/viscoelastic2/Kernels/Time.cpp
Equations/elastic/Model/Setup.cpp
Equations/elastic/Kernels/Neighbor.cpp
Equations/elastic/Kernels/Local.cpp
Equations/elastic/Kernels/Time.cpp
Equations/elastic/Physics/InitialField.cpp
# Equations/viscoelastic/Model/Setup.cpp
# Equations/viscoelastic/Kernels/Neighbor.cpp
# Equations/viscoelastic/Kernels/Local.cpp
# Equations/viscoelastic/Kernels/Time.cpp
Solver/Simulator.cpp
Solver/FreeSurfaceIntegrator.cpp
Solver/Interoperability.cpp
Solver/time_stepping/MiniSeisSol.cpp
Solver/time_stepping/TimeCluster.cpp
Solver/time_stepping/TimeManager.cpp
Kernels/DynamicRupture.cpp
Kernels/Plasticity.cpp
Kernels/TimeCommon.cpp
Kernels/Receiver.cpp
SeisSol.cpp
SourceTerm/Manager.cpp

SourceTerm/NRFReader.cpp # if netCDF

SourceTerm/PointSource.cpp
Parallel/Pin.cpp
Parallel/MPI.cpp
Parallel/mpiC.cpp
Parallel/FaultMPI.cpp
Geometry/GambitReader.cpp

# TODO: Only if hdf5?
Geometry/PUMLReader.cpp

Geometry/MeshReaderFBinding.cpp
Geometry/MeshTools.cpp
Monitoring/FlopCounter.cpp
Monitoring/LoopStatistics.cpp
Reader/readparC.cpp
#Reader/StressReaderC.cpp
#Reader/AsagiModule.cpp
Checkpoint/Manager.cpp

# TODO: Only if mpi?
Checkpoint/mpio/Wavefield.cpp
Checkpoint/mpio/FaultAsync.cpp
Checkpoint/mpio/Fault.cpp
Checkpoint/mpio/WavefieldAsync.cpp

# Checkpoint/sionlib/Wavefield.cpp
# Checkpoint/sionlib/Fault.cpp

# TODO: Only if HDF5!
Checkpoint/h5/Wavefield.cpp
Checkpoint/h5/Fault.cpp

Checkpoint/Backend.cpp
Checkpoint/Fault.cpp
Checkpoint/posix/Wavefield.cpp
Checkpoint/posix/Fault.cpp
ResultWriter/AnalysisWriter.cpp
ResultWriter/FreeSurfaceWriterExecutor.cpp
ResultWriter/PostProcessor.cpp
ResultWriter/FaultWriterC.cpp
ResultWriter/ReceiverWriter.cpp
ResultWriter/FaultWriterExecutor.cpp
ResultWriter/FaultWriter.cpp
ResultWriter/WaveFieldWriter.cpp
ResultWriter/FreeSurfaceWriter.cpp
main.cpp

# Fortran:
Monitoring/bindMonitoring.f90
Geometry/mpiextractmesh.f90
Geometry/allocate_mesh.f90
Geometry/MeshReaderCBinding.f90
Solver/close_seissol.f90
Solver/calc_deltat.f90
Solver/mpiexchangevalues.f90
Solver/prak_clif_mod.f90
Solver/calc_seissol.f90
Solver/f_ctof_bind_interoperability.f90
Solver/plasticity.f90
Solver/f_ftoc_bind_interoperability.f90
Numerical_aux/quadpoints.f90
Numerical_aux/jacobinormal.f90
Numerical_aux/convertxieta2xy.f90
Numerical_aux/create_fault_rotationmatrix.f90
Numerical_aux/trilinearinterpolation.f90
Numerical_aux/typesdef.f90
Numerical_aux/dgbasis.f90
Numerical_aux/gauss.f90
Numerical_aux/operators.f90
Modules/ModulesF.f90
seissolxx.f90
Physics/ini_model.f90
Physics/Evaluate_friction_law.f90
Physics/ini_model_DR.f90
# Reader/StressReaderF.f90
Reader/readpar.f90
Reader/read_backgroundstress.f90
ResultWriter/inioutput_seissol.f90
ResultWriter/magnitude_output.f90
# ResultWriter/receiver_hdf.f90
ResultWriter/output_rupturefront.f90
ResultWriter/ini_faultoutput.f90
ResultWriter/energies.f90
ResultWriter/FaultWriterF.f90
ResultWriter/faultoutput.f90
ResultWriter/common_fault_receiver.f90
ResultWriter/receiver.f90
Initializer/dg_setup.f90
Initializer/ini_optionalfields.f90
Initializer/ini_seissol.f90
Parallel/mpiF.f90
)


# Libs
find_package(yaml-cpp REQUIRED)

# TODO: Write export for impalaJIT
target_link_libraries(SeisSol PUBLIC yaml-cpp
			      ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/ImpalaJIT/libimpalajit.a 
			      )

# todo: if parallelization == omp/hybrid
include(FindOpenMP)
target_link_libraries(SeisSol PUBLIC OpenMP::OpenMP_CXX)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

include(FindMPI) # TODO: Doesnt work yet
target_include_directories(SeisSol PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(SeisSol PUBLIC ${MPI_CXX_LIBRARIES})
target_link_libraries(SeisSol PUBLIC ${MPI_C_LIBRARIES})
target_include_directories(SeisSol PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(SeisSol PUBLIC ${MPI_Fortran_LIBRARIES})

# TODO: HDF5 should be optional?
# TODO: Check constraints, seee
# https://github.com/Kitware/CMake/blob/master/Modules/FindHDF5.cmake
find_package(HDF5 REQUIRED)
target_include_directories(SeisSol PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(SeisSol PUBLIC ${HDF5_LIBRARIES})

find_package(NetCDF REQUIRED)
target_include_directories(SeisSol PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(SeisSol PUBLIC NetCDF::NetCDF)

# Parmetis
find_package(ParMETIS REQUIRED)
target_include_directories(SeisSol PUBLIC ${PARMETIS_INCLUDE_DIRS})
target_link_libraries(SeisSol PUBLIC ${PARMETIS_LIBRARIES})

target_include_directories(SeisSol PUBLIC
				   ${CMAKE_CURRENT_SOURCE_DIR}
				   ${CMAKE_CURRENT_SOURCE_DIR}/../submodules
				   ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/async
				   # todo: use reexport from easi?
				   ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/easi/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/yaml-cpp/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/ImpalaJIT/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/yateto/include
				   ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/glm
				   # todo: generate code in correct folder
				   ${CMAKE_CURRENT_SOURCE_DIR}/generated_code
				   # todo: Switch kernels?
				   ${CMAKE_CURRENT_SOURCE_DIR}/Equations/elastic
				 
)

target_compile_definitions(SeisSol PRIVATE
				   REAL_SIZE=8
				   CONVERGENCE_ORDER=3
				   ALIGNMENT=64
				   NUMBER_OF_QUANTITIES=9
				   NUMBER_OF_RELAXATION_MECHANISMS=0
				   # TODO: Remove next one
				   VERSION_STRING=""
				   __SSE3__
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
				   OMP # todo
				   USE_MPI # todo
				   PARALLEL # same as USE_MPI
				   USE_HDF # todo
				   USE_METIS # todo
				   USE_NETCDF # todo
				   )



# TODO: MPI/OpenMP, metis, hdf5, netcdf, puml
