# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

find_package(Doxygen REQUIRED)

# somehow, the CMake Doxygen module does not go well with a pre-defined Doxyfile
# i.e. it'll generate the data in the wrong place then
# thus, we'll need to configure manually

configure_file(Doxyfile.in Doxyfile)

doxygen_add_docs(seissol-doxygen
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/app
    ALL
    CONFIG_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
