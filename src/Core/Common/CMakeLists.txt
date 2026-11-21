# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-common OBJECT
    ConfigRegistry.cpp
    ConfigValue.cpp
    Environment.cpp
    Filesystem.cpp
    IntegerMaskParser.cpp
    Json.cpp
    Logging.cpp
    RunConfigs.cpp
    RunId.cpp
)

# The numbering of the types (TypeId.h) must exist once in the process. The library of the GPU
# code looks up variables of the storages too, but does not link the rest of SeisSol: in a GPU
# build, the numbering lives in that library, and the rest of SeisSol uses it from there.
if (NOT WITH_GPU)
  target_sources(seissol-common PRIVATE TypeId.cpp)
endif()

target_link_libraries(seissol-common PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-common)
