# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-plasticity OBJECT
    Setup/PlasticityPointGenerator.cpp
    Solver/Plasticity.cpp
    Solver/PlasticityRecorder.cpp
    Solver/PlasticityFeature.cpp
    Solver/PlasticityStage.cpp
)

target_link_libraries(seissol-plasticity PUBLIC seissol-common-properties)
# a feature uses the wave propagation, and the core below it
target_link_libraries(seissol-plasticity PUBLIC seissol-wavepropagation)
target_link_libraries(seissol-lib PUBLIC seissol-plasticity)
