# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-initializer OBJECT
    CellLocalMatrices.cpp
    InitialFieldProjection.cpp
    MemoryManager.cpp
    ParameterDB.cpp
    PointMapper.cpp

    BatchRecorders/DynamicRuptureRecorder.cpp
    BatchRecorders/LocalIntegrationRecorder.cpp
    BatchRecorders/NeighIntegrationRecorder.cpp
    BatchRecorders/PlasticityRecorder.cpp

    InitProcedure/Init.cpp
    InitProcedure/InitIO.cpp
    InitProcedure/InitLayout.cpp
    InitProcedure/InitMesh.cpp
    InitProcedure/InitModel.cpp
    InitProcedure/InitSideConditions.cpp
    InitProcedure/Internal/Buckets.cpp
    InitProcedure/Internal/LtsSetup.cpp
    InitProcedure/Internal/MeshLayout.cpp

    Parameters/CubeGeneratorParameters.cpp
    Parameters/DRParameters.cpp
    Parameters/InitializationParameters.cpp
    Parameters/LtsParameters.cpp
    Parameters/MeshParameters.cpp
    Parameters/ModelParameters.cpp
    Parameters/OutputParameters.cpp
    Parameters/ParameterReader.cpp
    Parameters/SeisSolParameters.cpp
    Parameters/SourceParameters.cpp

    TimeStepping/ClusterLayout.cpp
    TimeStepping/GlobalTimestep.cpp
    TimeStepping/Halo.cpp
    TimeStepping/LtsWeights/LtsWeights.cpp
    TimeStepping/LtsWeights/WeightsModels.cpp
)

target_link_libraries(seissol-initializer PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-initializer)
