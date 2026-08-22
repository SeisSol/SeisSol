# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-initializer OBJECT
    FaceMap.cpp
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
    InitProcedure/Internal/Boundary.cpp
    InitProcedure/Internal/Buckets.cpp
    InitProcedure/Internal/LtsSetup.cpp
    InitProcedure/Internal/MeshLayout.cpp
    InitProcedure/Internal/Recording.cpp
    InitProcedure/Internal/Scratchpads.cpp


    Model/BoundaryMappings.cpp
    Model/CellLocalMatrices.cpp
    Model/DynamicRuptureMatrices.cpp

    Parameters/CubeGeneratorParameters.cpp
    Parameters/DatafieldParameters.cpp
    Parameters/DRParameters.cpp
    Parameters/InitializationParameters.cpp
    Parameters/LtsParameters.cpp
    Parameters/MeshParameters.cpp
    Parameters/ModelParameters.cpp
    Parameters/OutputParameters.cpp
    Parameters/ParameterReader.cpp
    Parameters/SeisSolParameters.cpp
    Parameters/SourceParameters.cpp

    Clustering/ClusterHistogram.cpp
    Clustering/ClusterLadder.cpp
    Clustering/ClusterSmoother.cpp
    Clustering/Clustering.cpp
    Clustering/ClusteringCost.cpp
    Clustering/ClusteringEvaluator.cpp
    Clustering/GridLadderSearch.cpp
    Clustering/LadderOptimizer.cpp
    Clustering/LatticeDpSearch.cpp
    Clustering/TimestepHistogram.cpp
    Clustering/VertexWeights/VertexWeightModel.cpp
    Clustering/VertexWeights/WeightsModels.cpp

    TimeStepping/ClusterLayout.cpp
    TimeStepping/GlobalTimestep.cpp
    TimeStepping/Halo.cpp
)

target_link_libraries(seissol-initializer PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-initializer)
