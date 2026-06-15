# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-reader OBJECT
    AsagiModule.cpp
    AsagiReader.cpp
)

target_link_libraries(seissol-reader PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-reader)
