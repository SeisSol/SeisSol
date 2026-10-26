# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Smoke tests: end-to-end checks that run the real binaries rather than the
# doctest modules. They are labelled `smoke` and can be selected or skipped
# with `ctest -L smoke` / `ctest -LE smoke`.
#
# The checking itself lives in scripts/validate/smoke.py. Python3 is already a
# hard requirement of the build (see the top-level CMakeLists.txt), so this
# adds no dependency, and it keeps assertions about the binaries' output in a
# place where they can be unit-tested and linted like the rest of the tooling.
#
# TESTING_COMMAND is prepended to the command, matching what the doctest
# modules get through their CROSSCOMPILING_EMULATOR property. That property
# does not apply here, because the test command is the driver and not the
# target itself.

get_filename_component(SEISSOL_SMOKE_DRIVER
    "${CMAKE_CURRENT_LIST_DIR}/../scripts/validate/smoke.py" ABSOLUTE)

if (NOT EXISTS "${SEISSOL_SMOKE_DRIVER}")
    message(FATAL_ERROR "smoke-test driver not found: ${SEISSOL_SMOKE_DRIVER}")
endif()

# seissol_add_smoke_test(NAME <name> TARGET <target>
#                        [ARGS <arg>...]
#                        [EXPECT_FAILURE]
#                        [EXPECT_OUTPUT <regex>]
#                        [TIMEOUT <seconds>]
#                        [LABELS <label>...])
#
# Runs the target and asserts how it terminated.
function(seissol_add_smoke_test)
    cmake_parse_arguments(SMOKE
        "EXPECT_FAILURE"
        "NAME;TARGET;EXPECT_OUTPUT;TIMEOUT"
        "ARGS;LABELS"
        ${ARGN})

    if (NOT SMOKE_NAME OR NOT SMOKE_TARGET)
        message(FATAL_ERROR "seissol_add_smoke_test: NAME and TARGET are required")
    endif()

    set(_driver_args "")
    if (SMOKE_EXPECT_FAILURE)
        list(APPEND _driver_args --expect-failure)
    endif()
    if (SMOKE_EXPECT_OUTPUT)
        list(APPEND _driver_args --expect-output "${SMOKE_EXPECT_OUTPUT}")
    endif()
    if (SMOKE_TIMEOUT)
        list(APPEND _driver_args --timeout "${SMOKE_TIMEOUT}")
    endif()

    separate_arguments(_launcher NATIVE_COMMAND "${TESTING_COMMAND}")

    add_test(NAME "${SMOKE_NAME}"
        COMMAND ${Python3_EXECUTABLE} "${SEISSOL_SMOKE_DRIVER}" run
                ${_driver_args}
                -- ${_launcher} "$<TARGET_FILE:${SMOKE_TARGET}>" ${SMOKE_ARGS})

    set_tests_properties("${SMOKE_NAME}" PROPERTIES LABELS "smoke;${SMOKE_LABELS}")
endfunction()

# seissol_add_proxy_smoke_test(NAME <name> TARGET <target> KERNEL <kernel>
#                              [CELLS <n>] [TIMESTEPS <n>]
#                              [TIMEOUT <seconds>]
#                              [LABELS <label>...])
#
# Runs the proxy in JSON output mode and checks the document for formal
# correctness: it has to parse, carry exactly the documented field set, and
# hold finite numbers in a plausible range. No performance assertions, so this
# holds on any machine.
function(seissol_add_proxy_smoke_test)
    cmake_parse_arguments(PROXY
        ""
        "NAME;TARGET;KERNEL;CELLS;TIMESTEPS;TIMEOUT"
        "LABELS"
        ${ARGN})

    if (NOT PROXY_NAME OR NOT PROXY_TARGET OR NOT PROXY_KERNEL)
        message(FATAL_ERROR
            "seissol_add_proxy_smoke_test: NAME, TARGET and KERNEL are required")
    endif()

    if (NOT PROXY_CELLS)
        set(PROXY_CELLS 10)
    endif()
    if (NOT PROXY_TIMESTEPS)
        set(PROXY_TIMESTEPS 1)
    endif()

    set(_driver_args --kernel "${PROXY_KERNEL}"
                     --cells "${PROXY_CELLS}"
                     --timesteps "${PROXY_TIMESTEPS}")
    if (PROXY_TIMEOUT)
        list(APPEND _driver_args --timeout "${PROXY_TIMEOUT}")
    endif()

    separate_arguments(_launcher NATIVE_COMMAND "${TESTING_COMMAND}")

    add_test(NAME "${PROXY_NAME}"
        COMMAND ${Python3_EXECUTABLE} "${SEISSOL_SMOKE_DRIVER}" proxy
                ${_driver_args}
                -- ${_launcher} "$<TARGET_FILE:${PROXY_TARGET}>")

    set_tests_properties("${PROXY_NAME}" PROPERTIES LABELS "smoke;${PROXY_LABELS}")
endfunction()
