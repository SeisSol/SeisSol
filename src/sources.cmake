# Source code
add_library(SeisSol-common-lib

${CMAKE_SOURCE_DIR}/src/Initializer/CellLocalMatrices.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/GlobalData.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/AbstractGhostTimeCluster.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/AbstractTimeCluster.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/ActorState.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/CommunicationManager.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/DirectGhostTimeCluster.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/GhostTimeClusterWithCopy.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/MiniSeisSol.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/TimeCluster.cpp
${CMAKE_SOURCE_DIR}/src/Solver/time_stepping/TimeManager.cpp

${CMAKE_SOURCE_DIR}/src/Kernels/DynamicRupture.cpp
${CMAKE_SOURCE_DIR}/src/Kernels/Plasticity.cpp
${CMAKE_SOURCE_DIR}/src/Kernels/TimeCommon.cpp
${CMAKE_SOURCE_DIR}/src/Kernels/Touch.cpp
${CMAKE_SOURCE_DIR}/src/Kernels/PointSourceClusterOnHost.cpp

${CMAKE_SOURCE_DIR}/src/Common/IntegerMaskParser.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/FrictionSolver.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/LinearSlipWeakening.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/NoFault.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/SourceTimeFunction.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/ThermalPressurization/ThermalPressurization.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Initializer/BaseDRInitializer.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Initializer/ImposedSlipRatesInitializer.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Initializer/LinearSlipWeakeningInitializer.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Initializer/RateAndStateInitializer.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Misc.cpp

${CMAKE_SOURCE_DIR}/src/Equations/elastic/Kernels/GravitationalFreeSurfaceBC.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/PointMapper.cpp
${CMAKE_SOURCE_DIR}/src/Modules/Module.cpp
${CMAKE_SOURCE_DIR}/src/Modules/Modules.cpp

${CMAKE_SOURCE_DIR}/src/Monitoring/FlopCounter.cpp
${CMAKE_SOURCE_DIR}/src/Monitoring/ActorStateStatistics.cpp
${CMAKE_SOURCE_DIR}/src/Monitoring/LoopStatistics.cpp
${CMAKE_SOURCE_DIR}/src/Monitoring/Stopwatch.cpp
${CMAKE_SOURCE_DIR}/src/Monitoring/Unit.cpp

${CMAKE_SOURCE_DIR}/src/Kernels/Receiver.cpp
${CMAKE_SOURCE_DIR}/src/Model/Common.cpp
${CMAKE_SOURCE_DIR}/src/Numerical/Functions.cpp
${CMAKE_SOURCE_DIR}/src/Numerical/Statistics.cpp
${CMAKE_SOURCE_DIR}/src/Parallel/Pin.cpp
${CMAKE_SOURCE_DIR}/src/Physics/InstantaneousTimeMirrorManager.cpp
${CMAKE_SOURCE_DIR}/src/Solver/Pipeline/DrTuner.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/ClusteringWriter.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/AsyncIO.cpp

${CMAKE_SOURCE_DIR}/src/SourceTerm/FSRMReader.cpp
${CMAKE_SOURCE_DIR}/src/SourceTerm/PointSource.cpp
${CMAKE_SOURCE_DIR}/src/SourceTerm/Manager.cpp

${CMAKE_SOURCE_DIR}/src/Solver/Simulator.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/AnalysisWriter.cpp
)

target_compile_options(SeisSol-common-lib PRIVATE -fPIC)

if (SHARED)
add_library(SeisSol-lib SHARED)
else()
add_library(SeisSol-lib STATIC)
endif()

target_sources(SeisSol-lib PRIVATE
${CMAKE_SOURCE_DIR}/src/ResultWriter/EnergyOutput.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/FreeSurfaceWriter.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/FreeSurfaceWriterExecutor.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/MiniSeisSolWriter.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/PostProcessor.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/ReceiverWriter.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/ThreadsPinningWriter.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/WaveFieldWriter.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/FaultWriter.cpp
${CMAKE_SOURCE_DIR}/src/ResultWriter/FaultWriterExecutor.cpp

${CMAKE_SOURCE_DIR}/src/DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Output/FaultRefiner/FaultRefiners.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Output/OutputAux.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Output/OutputManager.cpp
${CMAKE_SOURCE_DIR}/src/DynamicRupture/Output/ReceiverBasedOutput.cpp

${CMAKE_SOURCE_DIR}/src/Equations/poroelastic/Model/Datastructures.cpp

${CMAKE_SOURCE_DIR}/src/Geometry/MeshReader.cpp
${CMAKE_SOURCE_DIR}/src/Geometry/MeshTools.cpp

${CMAKE_SOURCE_DIR}/src/Initializer/InitProcedure/Init.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/InitProcedure/InitIO.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/InitProcedure/InitMesh.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/InitProcedure/InitModel.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/InitProcedure/InitIO.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/InitProcedure/InitSideConditions.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/InitialFieldProjection.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/InternalState.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/MemoryAllocator.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/MemoryManager.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/ParameterDB.cpp

${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/CubeGeneratorParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/DRParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/InitializationParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/LtsParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/MeshParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/ModelParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/OutputParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/ParameterReader.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/SeisSolParameters.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/Parameters/SourceParameters.cpp

${CMAKE_SOURCE_DIR}/src/Initializer/TimeStepping/GlobalTimestep.cpp
${CMAKE_SOURCE_DIR}/src/Initializer/TimeStepping/LtsLayout.cpp

${CMAKE_SOURCE_DIR}/src/Initializer/Tree/Lut.cpp

${CMAKE_SOURCE_DIR}/src/Numerical/ODEInt.cpp
${CMAKE_SOURCE_DIR}/src/Numerical/ODEVector.cpp
${CMAKE_SOURCE_DIR}/src/Numerical/Transformation.cpp

${CMAKE_SOURCE_DIR}/src/Physics/InitialField.cpp

${CMAKE_SOURCE_DIR}/src/SeisSol.cpp

${CMAKE_SOURCE_DIR}/src/Solver/FreeSurfaceIntegrator.cpp

${CMAKE_SOURCE_DIR}/src/Reader/AsagiModule.cpp
${CMAKE_SOURCE_DIR}/src/Reader/AsagiReader.cpp

${CMAKE_SOURCE_DIR}/src/Parallel/Runtime/StreamOMP.cpp
)

set(SYCL_DEPENDENT_SRC_FILES
  ${CMAKE_SOURCE_DIR}/src/DynamicRupture/Factory.cpp
  ${CMAKE_SOURCE_DIR}/src/Parallel/MPI.cpp
  PARENT_SCOPE
)

set(SYCL_ONLY_SRC_FILES
  ${CMAKE_SOURCE_DIR}/src/Parallel/Runtime/StreamSycl.cpp
  ${CMAKE_SOURCE_DIR}/src/Parallel/AcceleratorDevice.cpp
  ${CMAKE_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.cpp
  ${CMAKE_SOURCE_DIR}/src/Kernels/PointSourceClusterOnDevice.cpp
  PARENT_SCOPE
)

target_compile_options(SeisSol-common-properties INTERFACE ${EXTRA_CXX_FLAGS})

if (HDF5 AND MPI)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Geometry/PartitioningLib.cpp
    ${CMAKE_SOURCE_DIR}/src/Geometry/PUMLReader.cpp
    ${CMAKE_SOURCE_DIR}/src/Initializer/TimeStepping/LtsWeights/LtsWeights.cpp
    ${CMAKE_SOURCE_DIR}/src/Initializer/TimeStepping/LtsWeights/WeightsModels.cpp
    )
endif()

if (NETCDF)
  target_sources(SeisSol-common-lib PRIVATE ${CMAKE_SOURCE_DIR}/src/SourceTerm/NRFReader.cpp)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Geometry/NetcdfReader.cpp
    ${CMAKE_SOURCE_DIR}/src/Geometry/CubeGenerator.cpp
    )
endif()


# Eqations have to be set at compile time currently.
if ("${EQUATIONS}" STREQUAL "elastic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Equations/elastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/elastic/Kernels/Local.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/elastic/Kernels/Neighbor.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/elastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_SOURCE_DIR}/src/Equations/elastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_ELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Local.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Neighbor.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_VISCOELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic2")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Neighbor.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Local.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_SOURCE_DIR}/src/Equations/viscoelastic2)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_VISCOELASTIC2)

elseif ("${EQUATIONS}" STREQUAL "anisotropic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Equations/anisotropic/Kernels/DirichletBoundary.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Neighbor.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Local.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_SOURCE_DIR}/src/Equations/anisotropic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_ANISOTROPIC)

elseif ("${EQUATIONS}" STREQUAL "poroelastic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Neighbor.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Local.cpp
    ${CMAKE_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_SOURCE_DIR}/src/Equations/poroelastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_STP)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_POROELASTIC)
endif()

target_include_directories(SeisSol-common-properties INTERFACE
        ${CMAKE_SOURCE_DIR}/src/Initializer/BatchRecorders)


# GPU code
if (WITH_GPU)
  target_sources(SeisSol-lib PRIVATE
          ${CMAKE_SOURCE_DIR}/src/Initializer/BatchRecorders/LocalIntegrationRecorder.cpp
          ${CMAKE_SOURCE_DIR}/src/Initializer/BatchRecorders/NeighIntegrationRecorder.cpp
          ${CMAKE_SOURCE_DIR}/src/Initializer/BatchRecorders/PlasticityRecorder.cpp
          ${CMAKE_SOURCE_DIR}/src/Initializer/BatchRecorders/DynamicRuptureRecorder.cpp)

  set(SEISSOL_DEVICE_INCLUDE ${DEVICE_INCLUDE_DIRS}
                             ${CMAKE_SOURCE_DIR}/submodules/yateto/include
                             ${CMAKE_BINARY_DIR}/generated-code
                             ${CMAKE_BINARY_DIR}/src
                             ${CMAKE_SOURCE_DIR}/src)

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
else()

endif()

