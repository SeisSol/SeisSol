# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Source code
add_library(seissol-kernel-lib

# run YATeTo first, since kernel.cpp usually takes really long

# kernel.cpp usually takes the longest
# (for CPUs, at least; for GPUs, we have a different library alltogether)
${CMAKE_BINARY_DIR}/codegen/GeneratedCode/kernel.cpp
${CMAKE_BINARY_DIR}/codegen/GeneratedCode/tensor.cpp
${CMAKE_BINARY_DIR}/codegen/GeneratedCode/subroutine.cpp
${CMAKE_BINARY_DIR}/codegen/GeneratedCode/init.cpp
)

# Generated code does only work without red-zone.
if (HAS_REDZONE)
  set_source_files_properties(
      ${CMAKE_BINARY_DIR}/codegen/GeneratedCode/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone
  )
endif()

target_compile_options(seissol-kernel-lib PRIVATE -fPIC)

if (SHARED)
  add_library(seissol-lib SHARED)
else()
  add_library(seissol-lib STATIC)
endif()

target_sources(seissol-lib PRIVATE

Initializer/CellLocalMatrices.cpp
Memory/GlobalData.cpp
Solver/TimeStepping/AbstractGhostTimeCluster.cpp
Solver/TimeStepping/AbstractTimeCluster.cpp
Solver/TimeStepping/ActorState.cpp
Solver/TimeStepping/CommunicationManager.cpp
Solver/TimeStepping/DirectGhostTimeCluster.cpp
Solver/TimeStepping/GhostTimeClusterWithCopy.cpp
Solver/TimeStepping/TimeCluster.cpp
Solver/TimeStepping/TimeManager.cpp

Kernels/DynamicRupture.cpp
Kernels/Plasticity.cpp
Kernels/TimeCommon.cpp
Kernels/Touch.cpp
Kernels/PointSourceClusterOnHost.cpp

Common/Filesystem.cpp
Common/IntegerMaskParser.cpp
DynamicRupture/FrictionLaws/FrictionSolver.cpp
DynamicRupture/FrictionLaws/CpuImpl/LinearSlipWeakening.cpp
DynamicRupture/FrictionLaws/CpuImpl/NoFault.cpp
DynamicRupture/FrictionLaws/CpuImpl/SourceTimeFunction.cpp
DynamicRupture/FrictionLaws/CpuImpl/ThermalPressurization/ThermalPressurization.cpp
DynamicRupture/Initializer/BaseDRInitializer.cpp
DynamicRupture/Initializer/ImposedSlipRatesInitializer.cpp
DynamicRupture/Initializer/LinearSlipWeakeningInitializer.cpp
DynamicRupture/Initializer/RateAndStateInitializer.cpp
DynamicRupture/Misc.cpp

Equations/anisotropic/Model/Datastructures.cpp
Equations/poroelastic/Model/Datastructures.cpp
Kernels/LinearCK/GravitationalFreeSurfaceBC.cpp

Initializer/InitialFieldProjection.cpp
Initializer/PointMapper.cpp
Modules/Module.cpp
Modules/Modules.cpp

Monitoring/FlopCounter.cpp
Monitoring/ActorStateStatistics.cpp
Monitoring/LoopStatistics.cpp
Monitoring/Stopwatch.cpp
Monitoring/Unit.cpp

Kernels/Receiver.cpp
Model/Common.cpp
Numerical/Functions.cpp
Numerical/Statistics.cpp
Parallel/Pin.cpp
Parallel/HelperThread.cpp
Physics/InstantaneousTimeMirrorManager.cpp
ResultWriter/ClusteringWriter.cpp
ResultWriter/AsyncIO.cpp

SourceTerm/FSRMReader.cpp
SourceTerm/PointSource.cpp
SourceTerm/Manager.cpp

Solver/Simulator.cpp
ResultWriter/AnalysisWriter.cpp

Initializer/TimeStepping/ClusterLayout.cpp
Solver/TimeStepping/HaloCommunication.cpp
Initializer/TimeStepping/Halo.cpp

DynamicRupture/Factory.cpp
Parallel/MPI.cpp

Parallel/OpenMP.cpp
Parallel/Runtime/Stream.cpp
Parallel/DataCollector.cpp

Solver/Estimator.cpp

ResultWriter/EnergyOutput.cpp
ResultWriter/FreeSurfaceWriter.cpp
ResultWriter/FreeSurfaceWriterExecutor.cpp
ResultWriter/MiniSeisSolWriter.cpp
ResultWriter/PostProcessor.cpp
ResultWriter/ReceiverWriter.cpp
ResultWriter/ThreadsPinningWriter.cpp
ResultWriter/WaveFieldWriter.cpp
ResultWriter/FaultWriter.cpp
ResultWriter/FaultWriterExecutor.cpp

DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.cpp
DynamicRupture/Output/FaultRefiner/FaultRefiners.cpp
DynamicRupture/Output/OutputAux.cpp
DynamicRupture/Output/OutputManager.cpp
DynamicRupture/Output/ReceiverBasedOutput.cpp

Geometry/MeshReader.cpp
Geometry/MeshTools.cpp

Geometry/CellTransform.cpp

Initializer/InitProcedure/Init.cpp
Initializer/InitProcedure/InitIO.cpp
Initializer/InitProcedure/InitMesh.cpp
Initializer/InitProcedure/InitModel.cpp
Initializer/InitProcedure/InitIO.cpp
Initializer/InitProcedure/InitSideConditions.cpp
Memory/MemoryAllocator.cpp
Initializer/MemoryManager.cpp
Initializer/ParameterDB.cpp

Initializer/InitProcedure/InitLayout.cpp
Initializer/InitProcedure/Internal/Buckets.cpp
Initializer/InitProcedure/Internal/MeshLayout.cpp
Initializer/InitProcedure/Internal/LtsSetup.cpp

Initializer/Parameters/CubeGeneratorParameters.cpp
Initializer/Parameters/DRParameters.cpp
Initializer/Parameters/InitializationParameters.cpp
Initializer/Parameters/LtsParameters.cpp
Initializer/Parameters/MeshParameters.cpp
Initializer/Parameters/ModelParameters.cpp
Initializer/Parameters/OutputParameters.cpp
Initializer/Parameters/ParameterReader.cpp
Initializer/Parameters/SeisSolParameters.cpp
Initializer/Parameters/SourceParameters.cpp

Initializer/TimeStepping/GlobalTimestep.cpp

Numerical/ODEInt.cpp
Numerical/ODEVector.cpp
Numerical/Transformation.cpp

Physics/InitialField.cpp

SeisSol.cpp

Solver/FreeSurfaceIntegrator.cpp

Reader/AsagiModule.cpp
Reader/AsagiReader.cpp

Geometry/CubeGenerator.cpp
)

if (HDF5 AND MPI)
  target_sources(seissol-lib PRIVATE
    Geometry/PartitioningLib.cpp
    Geometry/PUMLReader.cpp
    Initializer/TimeStepping/LtsWeights/LtsWeights.cpp
    Initializer/TimeStepping/LtsWeights/WeightsModels.cpp
    )
endif()

if (NETCDF)
  target_sources(seissol-lib PRIVATE SourceTerm/NRFReader.cpp)
endif()


# Eqations have to be set at compile time currently.
if ("${EQUATIONS}" STREQUAL "elastic")
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Time.cpp
    )
  target_include_directories(seissol-common-properties INTERFACE Equations/elastic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_ELASTIC)

elseif ("${EQUATIONS}" STREQUAL "acoustic")
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Time.cpp
    )
  target_include_directories(seissol-common-properties INTERFACE Equations/acoustic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_ACOUSTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic")
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Time.cpp
    )
  target_include_directories(seissol-common-properties INTERFACE Equations/viscoelastic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_VISCOELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic2")
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCKAnelastic/Neighbor.cpp
    Kernels/LinearCKAnelastic/Local.cpp
    Kernels/LinearCKAnelastic/Time.cpp
  )
  target_include_directories(seissol-common-properties INTERFACE Equations/viscoelastic2)
  target_compile_definitions(seissol-common-properties INTERFACE USE_VISCOELASTIC2)

elseif ("${EQUATIONS}" STREQUAL "anisotropic")
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Time.cpp
  )
  target_include_directories(seissol-common-properties INTERFACE Equations/anisotropic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_ANISOTROPIC)

elseif ("${EQUATIONS}" STREQUAL "poroelastic")
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Time.cpp
    Kernels/STP/Time.cpp
  )
  target_include_directories(seissol-common-properties INTERFACE Equations/poroelastic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_POROELASTIC)
endif()


# GPU code
if (WITH_GPU)
  target_sources(seissol-lib PRIVATE
          Initializer/BatchRecorders/LocalIntegrationRecorder.cpp
          Initializer/BatchRecorders/NeighIntegrationRecorder.cpp
          Initializer/BatchRecorders/PlasticityRecorder.cpp
          Initializer/BatchRecorders/DynamicRuptureRecorder.cpp
          Parallel/AcceleratorDevice.cpp
          DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.cpp
          Kernels/PointSourceClusterOnDevice.cpp)

  # include cmake files will define seissol-device-lib target
  if ("${DEVICE_BACKEND}" STREQUAL "cuda" OR "${DEVICE_BACKEND}" STREQUAL "hip")
    set(DEVICE_SRC ${DEVICE_SRC}
      ${CMAKE_BINARY_DIR}/codegen/GeneratedCode/gpulike_subroutine.cpp
      Kernels/DeviceAux/cudahip/PlasticityAux.cpp
      Kernels/LinearCK/DeviceAux/cudahip/KernelsAux.cpp
      DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverCudaHip.cpp
      Kernels/PointSourceClusterCudaHip.cpp)
  elseif ("${DEVICE_BACKEND}" STREQUAL "hipsycl" OR "${DEVICE_BACKEND}" STREQUAL "acpp" OR "${DEVICE_BACKEND}" STREQUAL "oneapi")
    set(DEVICE_SRC ${DEVICE_SRC}
          ${CMAKE_BINARY_DIR}/codegen/GeneratedCode/gpulike_subroutine.cpp
          Kernels/DeviceAux/sycl/PlasticityAux.cpp
          Kernels/LinearCK/DeviceAux/sycl/KernelsAux.cpp
          DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverSycl.cpp
          Kernels/PointSourceClusterSycl.cpp)
  endif()

  make_device_lib(seissol-device-lib "${DEVICE_SRC}")

  target_link_libraries(seissol-device-lib PRIVATE seissol-common-properties)

  target_compile_options(seissol-device-lib PRIVATE -fPIC)
  if ("${EQUATIONS}" STREQUAL "elastic")
    target_compile_definitions(seissol-device-lib PRIVATE USE_ELASTIC)
  elseif ("${EQUATIONS}" STREQUAL "acoustic")
    target_compile_definitions(seissol-device-lib PRIVATE USE_ACOUSTIC)
  elseif ("${EQUATIONS}" STREQUAL "viscoelastic")
    target_compile_definitions(seissol-device-lib PRIVATE USE_VISCOELASTIC)
  elseif ("${EQUATIONS}" STREQUAL "viscoelastic2")
    target_compile_definitions(seissol-device-lib PRIVATE USE_VISCOELASTIC2)
  elseif ("${EQUATIONS}" STREQUAL "anisotropic")
    target_compile_definitions(seissol-device-lib PRIVATE USE_ANISOTROPIC)
  elseif ("${EQUATIONS}" STREQUAL "poroelastic")
    target_compile_definitions(seissol-device-lib PRIVATE USE_POROELASTIC)
  endif()
  target_include_directories(seissol-lib PRIVATE ${DEVICE_INCLUDE_DIRS})

  if ("${EQUATIONS}" STREQUAL "elastic")
    target_compile_definitions(seissol-device-lib PRIVATE USE_ELASTIC)
  endif()

  if (USE_DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
    target_compile_definitions(seissol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
  endif()
endif()
