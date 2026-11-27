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
    Solver/PlaneWaveOperator.cpp
    Solver/TimeCommon.cpp
    Solver/WaveCluster.cpp
    Solver/WaveClusterFactory.cpp
)

# the kernels of the solvers the configurations advance their cells with (cmake/configs.cmake);
# stp takes the local and neighbor kernels of linearck
if ("linearck" IN_LIST SEISSOL_SOLVERS_USED OR "stp" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-wavepropagation PRIVATE
    Solver/LinearCK/Local.cpp
    Solver/LinearCK/Neighbor.cpp
    )
endif()
if ("linearck" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-wavepropagation PRIVATE
    Solver/LinearCK/Time.cpp
    )
  target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_LINEARCK)
endif()
if ("linearckanelastic" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-wavepropagation PRIVATE
    Solver/LinearCKAnelastic/Neighbor.cpp
    Solver/LinearCKAnelastic/Local.cpp
    Solver/LinearCKAnelastic/Time.cpp
    )
  target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_LINEARCKANELASTIC)
endif()
if ("stp" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-wavepropagation PRIVATE
    Solver/STP/Time.cpp
    )
  target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_STP)
endif()

target_link_libraries(seissol-wavepropagation PUBLIC seissol-common-properties)
# the wave propagation uses the core
get_property(SEISSOL_CORE_TARGETS GLOBAL PROPERTY SEISSOL_CORE_TARGETS)
target_link_libraries(seissol-wavepropagation PUBLIC ${SEISSOL_CORE_TARGETS})
target_link_libraries(seissol-lib PUBLIC seissol-wavepropagation)
