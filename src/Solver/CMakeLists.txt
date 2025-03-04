# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-solver
    Simulator.cpp
    FreeSurfaceIntegrator.cpp
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
    Clustering/Communication/CommunicationFactory.cpp
    Clustering/Communication/CommunicationManager.cpp
    Clustering/Communication/DirectMPINeighborCluster.cpp
)

add_dependencies(seissol-solver seissol-codegen)
target_link_libraries(seissol-solver PUBLIC seissol-common-properties)
