# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

if (PROXY_PYBINDING)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        find_package(pybind11 2.6.2 REQUIRED)
    else()
        find_package(pybind11 2.3.0 REQUIRED)
    endif()

    pybind11_add_module(seissol_proxy_bindings Bindings.cpp)
    target_link_libraries(seissol_proxy_bindings PUBLIC seissol-lib)
    file(COPY proxy-runners DESTINATION ${CMAKE_BINARY_DIR})
endif()

# C/C++ proxy interface
add_executable(seissol-proxy-exe Main.cpp)
target_link_libraries(seissol-proxy-exe PRIVATE seissol-lib)
if (NOT NEW_BINARY_NAMING)
    set_target_properties(seissol-proxy-exe PROPERTIES OUTPUT_NAME "SeisSol_proxy_${EXE_NAME_PREFIX}")
else()
    set_target_properties(seissol-proxy-exe PROPERTIES OUTPUT_NAME "proxyseissol-${NAME_SUFFIX}")
endif()
install(TARGETS seissol-proxy-exe RUNTIME)

if (TESTING)
    # Each kernel is run through the JSON output mode, so the run itself and
    # the numbers it reports are checked in one go.
    foreach(kernel all all_dr ader localwoader local neigh neigh_dr godunov_dr)
        seissol_add_proxy_smoke_test(
            NAME "Proxy: ${kernel}"
            TARGET seissol-proxy-exe
            KERNEL ${kernel}
            LABELS proxy
        )
    endforeach()

    # The plain format is what a human reads; keep it from silently going away.
    seissol_add_smoke_test(
        NAME "Proxy: plain output"
        TARGET seissol-proxy-exe
        ARGS 10 1 all -f plain
        EXPECT_OUTPUT "PERFORMANCE SUMMARY"
        LABELS proxy
    )

    # Argument handling: both of these have to be rejected, not silently run
    # with a default.
    seissol_add_smoke_test(
        NAME "Proxy: no arguments"
        TARGET seissol-proxy-exe
        EXPECT_FAILURE
        LABELS proxy
    )
    seissol_add_smoke_test(
        NAME "Proxy: unknown kernel"
        TARGET seissol-proxy-exe
        ARGS 10 1 not-a-kernel
        EXPECT_FAILURE
        LABELS proxy
    )
endif()
