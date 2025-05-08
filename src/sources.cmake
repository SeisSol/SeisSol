# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-common-lib

Initializer/CellLocalMatrices.cpp
Memory/GlobalData.cpp
Solver/time_stepping/AbstractGhostTimeCluster.cpp
Solver/time_stepping/AbstractTimeCluster.cpp
Solver/time_stepping/ActorState.cpp
Solver/time_stepping/CommunicationManager.cpp
Solver/time_stepping/DirectGhostTimeCluster.cpp
Solver/time_stepping/GhostTimeClusterWithCopy.cpp
Solver/time_stepping/TimeCluster.cpp
Solver/time_stepping/TimeManager.cpp

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
Physics/InstantaneousTimeMirrorManager.cpp
ResultWriter/ClusteringWriter.cpp
ResultWriter/AsyncIO.cpp

SourceTerm/FSRMReader.cpp
SourceTerm/PointSource.cpp
SourceTerm/Manager.cpp

Solver/Simulator.cpp
ResultWriter/AnalysisWriter.cpp

DynamicRupture/Factory.cpp
Parallel/MPI.cpp
)

# target_link_options(seissol-common-lib PUBLIC seissol-kernel-lib)
target_compile_options(seissol-kernel-lib PRIVATE -fPIC)
target_compile_options(seissol-common-lib PRIVATE -fPIC)

if (SHARED)
add_library(seissol-lib SHARED)
else()
add_library(seissol-lib STATIC)
endif()

target_sources(seissol-lib PRIVATE

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

Initializer/InitProcedure/Init.cpp
Initializer/InitProcedure/InitIO.cpp
Initializer/InitProcedure/InitMesh.cpp
Initializer/InitProcedure/InitModel.cpp
Initializer/InitProcedure/InitIO.cpp
Initializer/InitProcedure/InitSideConditions.cpp
Initializer/InternalState.cpp
Memory/MemoryAllocator.cpp
Initializer/MemoryManager.cpp
Initializer/ParameterDB.cpp

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
Initializer/TimeStepping/LtsLayout.cpp

Memory/Tree/Lut.cpp

Numerical/ODEInt.cpp
Numerical/ODEVector.cpp
Numerical/Transformation.cpp

Physics/InitialField.cpp

SeisSol.cpp

Solver/FreeSurfaceIntegrator.cpp

Reader/AsagiModule.cpp
Reader/AsagiReader.cpp

Parallel/Runtime/StreamOMP.cpp
)

set(SYCL_ONLY_SRC_FILES
  Parallel/AcceleratorDevice.cpp
  DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.cpp
  Kernels/PointSourceClusterOnDevice.cpp)

target_compile_options(seissol-common-properties INTERFACE ${EXTRA_CXX_FLAGS})

if (HDF5 AND MPI)
  target_sources(seissol-lib PRIVATE
    Geometry/PartitioningLib.cpp
    Geometry/PUMLReader.cpp
    Initializer/TimeStepping/LtsWeights/LtsWeights.cpp
    Initializer/TimeStepping/LtsWeights/WeightsModels.cpp
    )
endif()

if (NETCDF)
  target_sources(seissol-common-lib PRIVATE SourceTerm/NRFReader.cpp)
  target_sources(seissol-lib PRIVATE
    Geometry/NetcdfReader.cpp
    Geometry/CubeGenerator.cpp
    )
endif()


# Eqations have to be set at compile time currently.
if ("${EQUATIONS}" STREQUAL "elastic")
  target_sources(seissol-common-lib PRIVATE
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Time.cpp
    )
  target_include_directories(seissol-common-properties INTERFACE Equations/elastic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_ELASTIC)

elseif ("${EQUATIONS}" STREQUAL "acoustic")
  target_sources(seissol-common-lib PRIVATE
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Time.cpp
    )
  target_include_directories(seissol-common-properties INTERFACE Equations/acoustic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_ACOUSTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic")
  target_sources(seissol-common-lib PRIVATE
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Time.cpp
    )
  target_include_directories(seissol-common-properties INTERFACE Equations/viscoelastic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_VISCOELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic2")
  target_sources(seissol-common-lib PRIVATE
    Kernels/LinearCKAnelastic/Neighbor.cpp
    Kernels/LinearCKAnelastic/Local.cpp
    Kernels/LinearCKAnelastic/Time.cpp
  )
  target_include_directories(seissol-common-properties INTERFACE Equations/viscoelastic2)
  target_compile_definitions(seissol-common-properties INTERFACE USE_VISCOELASTIC2)

elseif ("${EQUATIONS}" STREQUAL "anisotropic")
  target_sources(seissol-common-lib PRIVATE
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Time.cpp
  )
  target_include_directories(seissol-common-properties INTERFACE Equations/anisotropic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_ANISOTROPIC)

elseif ("${EQUATIONS}" STREQUAL "poroelastic")
  target_sources(seissol-common-lib PRIVATE
    Kernels/LinearCK/Neighbor.cpp
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Time.cpp
    Kernels/STP/Time.cpp
  )
  target_include_directories(seissol-common-properties INTERFACE Equations/poroelastic)
  target_compile_definitions(seissol-common-properties INTERFACE USE_STP)
  target_compile_definitions(seissol-common-properties INTERFACE USE_POROELASTIC)
endif()

target_include_directories(seissol-common-properties INTERFACE
        Initializer/BatchRecorders)


# GPU code
if (WITH_GPU)
  target_sources(seissol-lib PRIVATE
          Initializer/BatchRecorders/LocalIntegrationRecorder.cpp
          Initializer/BatchRecorders/NeighIntegrationRecorder.cpp
          Initializer/BatchRecorders/PlasticityRecorder.cpp
          Initializer/BatchRecorders/DynamicRuptureRecorder.cpp)


  set(SEISSOL_DEVICE_INCLUDE ${DEVICE_INCLUDE_DIRS}
                             ${CMAKE_SOURCE_DIR}/submodules/yateto/include
                             ${CMAKE_BINARY_DIR}/generated-code
                             ${CMAKE_BINARY_DIR}/src)

  # include cmake files will define seissol-device-lib target
  if ("${DEVICE_BACKEND}" STREQUAL "cuda")
    include(cuda.cmake)
  elseif ("${DEVICE_BACKEND}" STREQUAL "hip")
    include(hip.cmake)
  elseif ("${DEVICE_BACKEND}" STREQUAL "hipsycl" OR "${DEVICE_BACKEND}" STREQUAL "acpp" OR "${DEVICE_BACKEND}" STREQUAL "oneapi")
    include(sycl.cmake)
  endif()

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
    target_compile_definitions(seissol-device-lib PRIVATE USE_STP)
    target_compile_definitions(seissol-device-lib PRIVATE USE_POROELASTIC)
  endif()
  target_include_directories(seissol-lib PRIVATE ${DEVICE_INCLUDE_DIRS})

  if ("${EQUATIONS}" STREQUAL "elastic")
    target_compile_definitions(seissol-device-lib PRIVATE USE_ELASTIC)
  endif()
endif()
