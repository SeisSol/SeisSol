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

add_subdirectory(Common)

add_subdirectory(DynamicRupture)

add_subdirectory(Equations)

add_subdirectory(Geometry)

add_subdirectory(Initializer)

add_subdirectory(IO)

add_subdirectory(Kernels)

add_subdirectory(Memory)

add_subdirectory(Model)

add_subdirectory(Modules)

add_subdirectory(Monitoring)

add_subdirectory(Numerical)

add_subdirectory(Parallel)

add_subdirectory(Physics)

add_subdirectory(Proxy)

add_subdirectory(Reader)

add_subdirectory(ResultWriter)

add_subdirectory(Solver)

add_subdirectory(SourceTerm)
