# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

find_package(Sphinx REQUIRED)

sphinx_add_docs(
  seissol-docs
  BUILDER html
  SOURCE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
