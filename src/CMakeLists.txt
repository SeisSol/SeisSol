# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

include(configprint.cmake)

include(sources.cmake)


# include subdirectories
# sort alphabetically
# leave one line in between (to facilitate Git merges)

add_subdirectory(Core)

add_subdirectory(Features)

add_subdirectory(Proxy)

add_subdirectory(Simulation)

add_subdirectory(WavePropagation)
