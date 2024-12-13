# Source code
add_library(SeisSol-kernel-lib

# do YATeTo first, since kernel.cpp usually takes really long

# kernel.cpp usually takes the longest
# (for CPUs, at least; for GPUs, we have a different library alltogether)
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/kernel.cpp
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/tensor.cpp
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/subroutine.cpp
${CMAKE_CURRENT_BINARY_DIR}/src/generated_code/init.cpp
)

add_library(SeisSol-common-lib

src/Initializer/CellLocalMatrices.cpp
src/Initializer/GlobalData.cpp
src/Solver/time_stepping/AbstractGhostTimeCluster.cpp
src/Solver/time_stepping/AbstractTimeCluster.cpp
src/Solver/time_stepping/ActorState.cpp
src/Solver/time_stepping/CommunicationManager.cpp
src/Solver/time_stepping/DirectGhostTimeCluster.cpp
src/Solver/time_stepping/GhostTimeClusterWithCopy.cpp
src/Solver/time_stepping/MiniSeisSol.cpp
src/Solver/time_stepping/TimeCluster.cpp
src/Solver/time_stepping/TimeManager.cpp

src/Kernels/DynamicRupture.cpp
src/Kernels/Plasticity.cpp
src/Kernels/TimeCommon.cpp
src/Kernels/Touch.cpp
src/Kernels/PointSourceClusterOnHost.cpp

src/Common/Filesystem.cpp
src/Common/IntegerMaskParser.cpp
src/DynamicRupture/FrictionLaws/FrictionSolver.cpp
src/DynamicRupture/FrictionLaws/LinearSlipWeakening.cpp
src/DynamicRupture/FrictionLaws/NoFault.cpp
src/DynamicRupture/FrictionLaws/SourceTimeFunction.cpp
src/DynamicRupture/FrictionLaws/ThermalPressurization/ThermalPressurization.cpp
src/DynamicRupture/Initializer/BaseDRInitializer.cpp
src/DynamicRupture/Initializer/ImposedSlipRatesInitializer.cpp
src/DynamicRupture/Initializer/LinearSlipWeakeningInitializer.cpp
src/DynamicRupture/Initializer/RateAndStateInitializer.cpp
src/DynamicRupture/Misc.cpp

src/Equations/elastic/Kernels/GravitationalFreeSurfaceBC.cpp
src/Initializer/InitialFieldProjection.cpp
src/Initializer/PointMapper.cpp
src/Modules/Module.cpp
src/Modules/Modules.cpp

src/Monitoring/FlopCounter.cpp
src/Monitoring/ActorStateStatistics.cpp
src/Monitoring/LoopStatistics.cpp
src/Monitoring/Stopwatch.cpp
src/Monitoring/Unit.cpp

src/Kernels/Receiver.cpp
src/Model/Common.cpp
src/Numerical/Functions.cpp
src/Numerical/Statistics.cpp
src/Parallel/Pin.cpp
src/Physics/InstantaneousTimeMirrorManager.cpp
src/ResultWriter/ClusteringWriter.cpp
src/ResultWriter/AsyncIO.cpp

src/SourceTerm/FSRMReader.cpp
src/SourceTerm/PointSource.cpp
src/SourceTerm/Manager.cpp

src/Solver/Simulator.cpp
src/ResultWriter/AnalysisWriter.cpp
)

# target_link_options(SeisSol-common-lib PUBLIC SeisSol-kernel-lib)
target_compile_options(SeisSol-kernel-lib PRIVATE -fPIC)
target_compile_options(SeisSol-common-lib PRIVATE -fPIC)

if (SHARED)
add_library(SeisSol-lib SHARED)
else()
add_library(SeisSol-lib STATIC)
endif()

target_sources(SeisSol-lib PRIVATE
src/ResultWriter/EnergyOutput.cpp
src/ResultWriter/FreeSurfaceWriter.cpp
src/ResultWriter/FreeSurfaceWriterExecutor.cpp
src/ResultWriter/MiniSeisSolWriter.cpp
src/ResultWriter/PostProcessor.cpp
src/ResultWriter/ReceiverWriter.cpp
src/ResultWriter/ThreadsPinningWriter.cpp
src/ResultWriter/WaveFieldWriter.cpp
src/ResultWriter/FaultWriter.cpp
src/ResultWriter/FaultWriterExecutor.cpp

src/DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.cpp
src/DynamicRupture/Output/FaultRefiner/FaultRefiners.cpp
src/DynamicRupture/Output/OutputAux.cpp
src/DynamicRupture/Output/OutputManager.cpp
src/DynamicRupture/Output/ReceiverBasedOutput.cpp

src/Equations/poroelastic/Model/Datastructures.cpp

src/Geometry/MeshReader.cpp
src/Geometry/MeshTools.cpp

src/Initializer/InitProcedure/Init.cpp
src/Initializer/InitProcedure/InitIO.cpp
src/Initializer/InitProcedure/InitMesh.cpp
src/Initializer/InitProcedure/InitModel.cpp
src/Initializer/InitProcedure/InitIO.cpp
src/Initializer/InitProcedure/InitSideConditions.cpp
src/Initializer/InternalState.cpp
src/Initializer/MemoryAllocator.cpp
src/Initializer/MemoryManager.cpp
src/Initializer/ParameterDB.cpp

src/Initializer/Parameters/CubeGeneratorParameters.cpp
src/Initializer/Parameters/DRParameters.cpp
src/Initializer/Parameters/InitializationParameters.cpp
src/Initializer/Parameters/LtsParameters.cpp
src/Initializer/Parameters/MeshParameters.cpp
src/Initializer/Parameters/ModelParameters.cpp
src/Initializer/Parameters/OutputParameters.cpp
src/Initializer/Parameters/ParameterReader.cpp
src/Initializer/Parameters/SeisSolParameters.cpp
src/Initializer/Parameters/SourceParameters.cpp

src/Initializer/TimeStepping/GlobalTimestep.cpp
src/Initializer/TimeStepping/LtsLayout.cpp

src/Initializer/Tree/Lut.cpp

src/Numerical/ODEInt.cpp
src/Numerical/ODEVector.cpp
src/Numerical/Transformation.cpp

src/Physics/InitialField.cpp

src/SeisSol.cpp

src/Solver/FreeSurfaceIntegrator.cpp

src/Reader/AsagiModule.cpp
src/Reader/AsagiReader.cpp

src/Parallel/Runtime/StreamOMP.cpp
)

set(SYCL_DEPENDENT_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/Factory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Parallel/MPI.cpp
)

set(SYCL_ONLY_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Parallel/Runtime/StreamSycl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Parallel/AcceleratorDevice.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/PointSourceClusterOnDevice.cpp)

target_compile_options(SeisSol-common-properties INTERFACE ${EXTRA_CXX_FLAGS})
target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/src/generated_code)

if (HDF5 AND MPI)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/PartitioningLib.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/PUMLReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/TimeStepping/LtsWeights/LtsWeights.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/TimeStepping/LtsWeights/WeightsModels.cpp
    )
endif()

if (NETCDF)
  target_sources(SeisSol-common-lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceTerm/NRFReader.cpp)
  target_sources(SeisSol-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/NetcdfReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/CubeGenerator.cpp
    )
endif()


# Eqations have to be set at compile time currently.
if ("${EQUATIONS}" STREQUAL "elastic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_ELASTIC)

elseif ("${EQUATIONS}" STREQUAL "acoustic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/acoustic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/acoustic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/acoustic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/acoustic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/acoustic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_ACOUSTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_VISCOELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic2")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_VISCOELASTIC2)

elseif ("${EQUATIONS}" STREQUAL "anisotropic")
  target_sources(SeisSol-common-lib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-common-properties INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic)
  target_compile_definitions(SeisSol-common-properties INTERFACE USE_ANISOTROPIC)

elseif ("${EQUATIONS}" STREQUAL "poroelastic")
  target_sources(SeisSol-common-lib PRIVATE
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
                             ${CMAKE_BINARY_DIR}/src
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

add_subdirectory(src/IO)
target_link_libraries(SeisSol-lib PUBLIC seissol-io)
