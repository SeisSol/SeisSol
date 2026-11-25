# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-initialconditions OBJECT
    EasiFields.cpp
    Output/AnalysisWriter.cpp
    Physics/Scenario/Registry.cpp
    Physics/Scenario/Scenarios.cpp
    Setup/InitialFieldProjection.cpp
)

target_link_libraries(seissol-initialconditions PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-initialconditions)
