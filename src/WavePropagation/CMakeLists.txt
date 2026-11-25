# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-wavepropagation OBJECT
    Output/Checkpoint.cpp
    Output/DofSync.cpp
    Output/FreeSurfaceIntegrator.cpp
    Output/IntegrationFeature.cpp
    Output/IntegrationStage.cpp
    Output/Receiver.cpp
    Output/ReceiverWriter.cpp
    Physics/Common.cpp
    Physics/anisotropic/Datastructures.cpp
    Physics/poroelastic/Datastructures.cpp
    Setup/BoundaryMappings.cpp
    Setup/CellLocalMatrices.cpp
    Setup/DirichletCondition.cpp
    Setup/MaterialParameterDB.cpp
    Setup/Timesteps.cpp
    Solver/ConfigBoundary.cpp
    Solver/LocalIntegrationRecorder.cpp
    Solver/NeighIntegrationRecorder.cpp
    Solver/TimeCommon.cpp
    Solver/WaveCluster.cpp
    Solver/WaveClusterFactory.cpp
)

target_link_libraries(seissol-wavepropagation PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-wavepropagation)
