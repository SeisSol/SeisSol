# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_executable(seissol-bin Main.cpp)
if (NOT NEW_BINARY_NAMING)
  set_target_properties(seissol-bin PROPERTIES OUTPUT_NAME "SeisSol_${EXE_NAME_PREFIX}")
else()
  set_target_properties(seissol-bin PROPERTIES OUTPUT_NAME "seissol-${NAME_SUFFIX}")
endif()

target_link_libraries(seissol-bin PUBLIC seissol-lib)

install(TARGETS seissol-bin RUNTIME)

if (TESTING)
    # Started without a parameter file, SeisSol has to get as far as reading it
    # and then bail out through logError(). This covers static initialisation,
    # MPI startup and argument parsing for every equation/precision variant,
    # which is where a broken build usually shows first.
    seissol_add_smoke_test(
        NAME "Smoke: no parameter file"
        TARGET seissol-bin
        EXPECT_FAILURE
        EXPECT_OUTPUT "Error while reading the parameter file"
    )
endif()
