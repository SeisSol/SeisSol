# Source code
add_library(SeisSol-lib

# do YATeTo first, since kernel.cpp usually takes really long

# kernel.cpp usually takes the longest
# (for CPUs, at least; for GPUs, we have a different library alltogether)
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/kernel.cpp
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/tensor.cpp
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/subroutine.cpp
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/init.cpp

src/Initializer/ParameterDB.cpp
src/Initializer/PointMapper.cpp
src/Initializer/GlobalData.cpp
src/Initializer/InternalState.cpp
src/Initializer/MemoryAllocator.cpp
src/Initializer/CellLocalMatrices.cpp

src/Initializer/time_stepping/LtsLayout.cpp
src/Initializer/time_stepping/LtsParameters.cpp
src/Initializer/time_stepping/GlobalTimestep.cpp
src/Initializer/tree/Lut.cpp
src/Initializer/MemoryManager.cpp
src/Initializer/InitialFieldProjection.cpp
src/Initializer/InputParameters.cpp

src/Initializer/InitProcedure/InitMesh.cpp
src/Initializer/InitProcedure/InitModel.cpp
src/Initializer/InitProcedure/InitIO.cpp
src/Initializer/InitProcedure/InitSideConditions.cpp
src/Initializer/InitProcedure/Init.cpp

src/Modules/Modules.cpp
src/Model/common.cpp
src/Numerical_aux/Functions.cpp
src/Numerical_aux/Transformation.cpp
src/Numerical_aux/Statistics.cpp

src/Solver/Simulator.cpp
src/Solver/FreeSurfaceIntegrator.cpp

src/Solver/time_stepping/AbstractTimeCluster.cpp
src/Solver/time_stepping/ActorState.cpp
src/Solver/time_stepping/MiniSeisSol.cpp
src/Solver/time_stepping/TimeCluster.cpp
src/Solver/time_stepping/AbstractGhostTimeCluster.cpp
src/Solver/time_stepping/DirectGhostTimeCluster.cpp
src/Solver/time_stepping/GhostTimeClusterWithCopy.cpp
src/Solver/time_stepping/CommunicationManager.cpp

src/Solver/time_stepping/TimeManager.cpp
src/Solver/Pipeline/DrTuner.cpp
src/Kernels/DynamicRupture.cpp
src/Kernels/Plasticity.cpp
src/Kernels/TimeCommon.cpp
src/Kernels/Receiver.cpp
src/Kernels/Touch.cpp
src/SeisSol.cpp
src/Parallel/Pin.cpp

src/Geometry/MeshTools.cpp
src/Geometry/MeshReader.cpp
src/Monitoring/FlopCounter.cpp
src/Monitoring/LoopStatistics.cpp
src/Monitoring/ActorStateStatistics.cpp
src/Monitoring/Stopwatch.cpp
src/Monitoring/Unit.cpp

src/Checkpoint/Manager.cpp

src/Checkpoint/Backend.cpp
src/Checkpoint/Fault.cpp
src/Checkpoint/posix/Wavefield.cpp
src/Checkpoint/posix/Fault.cpp
src/ResultWriter/AnalysisWriter.cpp
src/ResultWriter/MiniSeisSolWriter.cpp
src/ResultWriter/ClusteringWriter.cpp
src/ResultWriter/EnergyOutput.cpp
src/ResultWriter/ThreadsPinningWriter.cpp
src/ResultWriter/FreeSurfaceWriterExecutor.cpp
src/ResultWriter/PostProcessor.cpp
src/ResultWriter/ReceiverWriter.cpp
src/ResultWriter/FaultWriterExecutor.cpp
src/ResultWriter/FaultWriter.cpp
src/ResultWriter/WaveFieldWriter.cpp
src/ResultWriter/FreeSurfaceWriter.cpp

src/Numerical_aux/ODEInt.cpp
src/Numerical_aux/ODEVector.cpp
src/Physics/Attenuation.cpp
src/Physics/InitialField.cpp
src/Physics/InstantaneousTimeMirrorManager.cpp

src/Equations/poroelastic/Model/datastructures.cpp
src/Equations/elastic/Kernels/GravitationalFreeSurfaceBC.cpp

src/Common/IntegerMaskParser.cpp


)

set(SYCL_DEPENDENT_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Model/common.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Parallel/MPI.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/FrictionSolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Misc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Factory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/SourceTimeFunction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/LinearSlipWeakening.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/NoFault.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/ThermalPressurization/ThermalPressurization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Initializers/BaseDRInitializer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Initializers/ImposedSlipRatesInitializer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Initializers/LinearSlipWeakeningInitializer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Initializers/RateAndStateInitializer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Output/OutputManager.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Output/ReceiverBasedOutput.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Output/FaultRefiner/FaultRefiners.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Output/OutputAux.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/PointSourceClusterOnHost.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceTerm/FSRMReader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceTerm/Manager.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceTerm/PointSource.cpp)

set(SYCL_ONLY_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Parallel/AcceleratorDevice.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/PointSourceClusterOnDevice.cpp)

target_compile_options(SeisSol-common-properties INTERFACE ${EXTRA_CXX_FLAGS})
target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code)

if (MPI)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/Wavefield.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/FaultAsync.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/Fault.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/WavefieldAsync.cpp
)
endif()

if (HDF5)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/h5/Wavefield.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/h5/Fault.cpp
    )
endif()

if (HDF5 AND MPI)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/PartitioningLib.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/PUMLReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/time_stepping/LtsWeights/LtsWeights.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/time_stepping/LtsWeights/WeightsModels.cpp
    )
endif()

if (NETCDF)
  list(APPEND SYCL_DEPENDENT_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceTerm/NRFReader.cpp)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/NetcdfReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/CubeGenerator.cpp
    )
endif()

if (ASAGI)
  target_sources(SeisSol-lib PRIVATE
    #todo:
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Reader/AsagiModule.cpp
    )
endif()


# Eqations have to be set at compile time currently.
if ("${EQUATIONS}" STREQUAL "elastic")
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_ELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic")
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_VISCOELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic2")
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_VISCOELASTIC2)

elseif ("${EQUATIONS}" STREQUAL "anisotropic")
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_ANISOTROPIC)

elseif ("${EQUATIONS}" STREQUAL "poroelastic")
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_STP)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_POROELASTIC)
endif()

target_include_directories(SeisSol-common-properties INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders)


# GPU code
if (WITH_GPU)
  target_sources(SeisSol-lib PRIVATE
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/LocalIntegrationRecorder.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/NeighIntegrationRecorder.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/PlasticityRecorder.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/DynamicRuptureRecorder.cpp)


  set(SEISSOL_DEVICE_INCLUDE ${DEVICE_INCLUDE_DIRS}
                             ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
                             ${CMAKE_BINARY_DIR}/src/generated_code
                             ${CMAKE_CURRENT_SOURCE_DIR}/src)

  # include cmake files will define SeisSol-device-lib target
  if ("${DEVICE_BACKEND}" STREQUAL "cuda")
    include(${CMAKE_SOURCE_DIR}/src/cuda.cmake)
  elseif ("${DEVICE_BACKEND}" STREQUAL "hip")
    include(${CMAKE_SOURCE_DIR}/src/hip.cmake)
  elseif ("${DEVICE_BACKEND}" STREQUAL "hipsycl" OR "${DEVICE_BACKEND}" STREQUAL "oneapi")
    include(${CMAKE_SOURCE_DIR}/src/sycl.cmake)
  endif()

  target_compile_options(SeisSol-device-lib PRIVATE -fPIC)
  target_include_directories(SeisSol-lib PRIVATE ${DEVICE_INCLUDE_DIRS})

  if ("${EQUATIONS}" STREQUAL "elastic")
    target_compile_definitions(SeisSol-device-lib PRIVATE USE_ELASTIC)
  endif()
endif()
