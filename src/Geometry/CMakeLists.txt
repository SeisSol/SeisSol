# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-geometry OBJECT
    CubeGenerator.cpp
    MeshReader.cpp
    MeshTools.cpp
    PartitioningLib.cpp
    PUMLReader.cpp
)

target_link_libraries(seissol-geometry PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-geometry)
