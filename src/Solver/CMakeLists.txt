# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

add_library(seissol-solver
    Simulator.cpp
    FreeSurfaceIntegrator.cpp
    Proxy/MiniSeisSol.cpp
    Clustering/TimeManager.cpp
    Clustering/ActorState.cpp
    Clustering/AbstractTimeCluster.cpp
    Clustering/Computation/TimeCluster.cpp
    Clustering/Computation/DynamicRuptureCluster.cpp
    Clustering/Computation/GhostCluster.cpp
    Clustering/Computation/CopyCluster.cpp
    Clustering/Communication/CCLNeighborCluster.cpp
    Clustering/Communication/CCLSetup.cpp
    Clustering/Communication/NeighborCluster.cpp
    Clustering/Communication/HaloHelper.cpp
    Clustering/Communication/CommunicationFactory.cpp
    Clustering/Communication/CommunicationManager.cpp
    Clustering/Communication/DirectMPINeighborCluster.cpp
)

add_dependencies(seissol-solver SeisSol-codegen)
target_link_libraries(seissol-solver PUBLIC SeisSol-common-properties)
