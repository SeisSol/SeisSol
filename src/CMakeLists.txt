# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

include(sources.cmake)

add_subdirectory(IO)
target_link_libraries(seissol-lib PUBLIC seissol-io)

add_subdirectory(Proxy)
