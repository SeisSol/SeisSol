# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-timestepping OBJECT
    Actor/AbstractTimeCluster.cpp
    Actor/ActorState.cpp
    Actor/StepParams.cpp
    Clustering/ClusterHistogram.cpp
    Clustering/ClusterLadder.cpp
    Clustering/ClusterSmoother.cpp
    Clustering/Clustering.cpp
    Clustering/ClusteringCost.cpp
    Clustering/ClusteringEvaluator.cpp
    Clustering/ClusteringWriter.cpp
    Clustering/GridLadderSearch.cpp
    Clustering/LadderOptimizer.cpp
    Clustering/LatticeDpSearch.cpp
    Clustering/TimestepHistogram.cpp
    Clustering/VertexWeights/MeasuredWeights.cpp
    Clustering/VertexWeights/VertexWeightModel.cpp
    Clustering/VertexWeights/WeightsModels.cpp
    Compute/CellCluster.cpp
    Compute/ClusterClock.cpp
    Compute/FaceCluster.cpp
    Halo/CommunicationManager.cpp
    Halo/GhostCluster.cpp
    Halo/HaloCheck.cpp
    Halo/Mpi/MpiHaloTransport.cpp
    Halo/Mpi/StagedMpiHaloTransport.cpp
    Halo/TaskHaloTransport.cpp
    HaloTransportFactory.cpp
    Layout/ClusterLayout.cpp
    Layout/GlobalTimestep.cpp
    Layout/Halo.cpp
    Plan/SuperStepRecorder.cpp
    Plan/TimeSteppingPlan.cpp
    TimeManager.cpp
)

target_link_libraries(seissol-timestepping PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-timestepping)
