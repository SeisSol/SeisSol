# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-solver OBJECT
    Estimator.cpp
    FreeSurfaceIntegrator.cpp
    Simulator.cpp
    TimeStepping/AbstractGhostTimeCluster.cpp
    TimeStepping/AbstractTimeCluster.cpp
    TimeStepping/ActorState.cpp
    TimeStepping/CommunicationManager.cpp
    TimeStepping/DirectGhostTimeCluster.cpp
    TimeStepping/GhostTimeClusterWithCopy.cpp
    TimeStepping/HaloCommunication.cpp
    TimeStepping/TimeCluster.cpp
    TimeStepping/TimeManager.cpp
)

target_link_libraries(seissol-solver PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-solver)
