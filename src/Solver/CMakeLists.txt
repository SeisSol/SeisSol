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
    TimeStepping/AbstractTimeCluster.cpp
    TimeStepping/ActorState.cpp
    TimeStepping/CclExchangeScheduler.cpp
    TimeStepping/CellCluster.cpp
    TimeStepping/ClusterClock.cpp
    TimeStepping/CommunicationManager.cpp
    TimeStepping/DynamicRuptureCluster.cpp
    TimeStepping/ExchangeScheduler.cpp
    TimeStepping/FaceCluster.cpp
    TimeStepping/GhostCluster.cpp
    TimeStepping/HaloTransport.cpp
    TimeStepping/MpiHaloTransport.cpp
    TimeStepping/StagedMpiHaloTransport.cpp
    TimeStepping/StepParams.cpp
    TimeStepping/SuperStepRecorder.cpp
    TimeStepping/TimeManager.cpp
    TimeStepping/TimeSteppingPlan.cpp
)

target_link_libraries(seissol-solver PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-solver)
