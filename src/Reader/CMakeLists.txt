# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-reader OBJECT
    Datafield/AsagiLite.cpp

    Datafield/AsagiModule.cpp
    Datafield/AsagiReader.cpp

    Scripting/EasiReader.cpp
    Scripting/LuaReader.cpp
    Scripting/ReaderBuilder.cpp
)

target_link_libraries(seissol-reader PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-reader)
