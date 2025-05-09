# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_executable(seissol-bin Main.cpp)
set_target_properties(seissol-bin PROPERTIES OUTPUT_NAME "SeisSol_${EXE_NAME_PREFIX}")

target_link_libraries(seissol-bin PUBLIC seissol-lib)

install(TARGETS seissol-bin RUNTIME)
