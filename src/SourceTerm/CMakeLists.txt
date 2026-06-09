# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-sourceterm OBJECT
    FSRMReader.cpp
    Manager.cpp
    PointSource.cpp
)

if (NETCDF)
    target_sources(seissol-sourceterm PRIVATE NRFReader.cpp)
endif()

target_link_libraries(seissol-sourceterm PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-sourceterm)
