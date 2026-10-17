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
    TimeStepping/Actor/AbstractTimeCluster.cpp
    TimeStepping/Actor/ActorState.cpp
    TimeStepping/Actor/StepParams.cpp
    TimeStepping/Compute/CellCluster.cpp
    TimeStepping/Compute/ClusterClock.cpp
    TimeStepping/Compute/DynamicRuptureCluster.cpp
    TimeStepping/Compute/FaceCluster.cpp
    TimeStepping/Halo/CommunicationManager.cpp
    TimeStepping/Halo/GhostCluster.cpp
    TimeStepping/Halo/Mpi/MpiHaloTransport.cpp
    TimeStepping/Halo/Mpi/StagedMpiHaloTransport.cpp
    TimeStepping/Halo/Stream/CclExchangeScheduler.cpp
    TimeStepping/Halo/Stream/ExchangeScheduler.cpp
    TimeStepping/Halo/Stream/ShmemExchangeScheduler.cpp
    TimeStepping/Halo/Stream/StreamExchangeScheduler.cpp
    TimeStepping/Halo/Stream/StreamMpiExchangeScheduler.cpp
    TimeStepping/HaloTransportFactory.cpp
    TimeStepping/Plan/SuperStepRecorder.cpp
    TimeStepping/Plan/TimeSteppingPlan.cpp
    TimeStepping/TimeManager.cpp
)

target_link_libraries(seissol-solver PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-solver)
