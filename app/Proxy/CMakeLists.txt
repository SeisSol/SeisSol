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
if (TARGET device-hsa-queue-hook)
  # the HSA queues of the stream runtime (SEISSOL_HSA_QUEUES, HIP only), as for seissol-bin
  target_link_libraries(seissol-proxy-exe PRIVATE device-hsa-queue-hook)
endif()
if (NOT NEW_BINARY_NAMING)
    set_target_properties(seissol-proxy-exe PROPERTIES OUTPUT_NAME "SeisSol_proxy_${EXE_NAME_PREFIX}")
else()
    set_target_properties(seissol-proxy-exe PROPERTIES OUTPUT_NAME "proxyseissol${NAME_SUFFIX}")
endif()
install(TARGETS seissol-proxy-exe RUNTIME)

if (TESTING)
    # Each kernel is run through the JSON output mode, so the run itself and
    # the numbers it reports are checked in one go: in the kernel mode (the
    # default), and the loops of the clusters in the cluster mode.
    foreach(kernel all all_dr local neigh neigh_dr godunov_dr ader localwoader plast plast_yield)
        seissol_add_proxy_smoke_test(
            NAME "Proxy: ${kernel}"
            TARGET seissol-proxy-exe
            KERNEL ${kernel}
            LABELS proxy
        )
    endforeach()
    foreach(kernel all all_dr local neigh neigh_dr godunov_dr)
        seissol_add_proxy_smoke_test(
            NAME "Proxy (clusters): ${kernel}"
            TARGET seissol-proxy-exe
            KERNEL ${kernel}
            OPTIONS --mode=cluster
            LABELS proxy
        )
    endforeach()
    foreach(plasticity elastic yield)
        seissol_add_proxy_smoke_test(
            NAME "Proxy (clusters, plasticity ${plasticity}): all"
            TARGET seissol-proxy-exe
            KERNEL all
            OPTIONS --mode=cluster --plasticity=${plasticity}
            LABELS proxy
        )
    endforeach()
    # the friction laws: alone in the kernel mode, with the faces in the cluster mode
    foreach(law 0 16 6 1058 3 4 7 103 103tp 33 34 35)
        seissol_add_proxy_smoke_test(
            NAME "Proxy (friction law ${law}): fl"
            TARGET seissol-proxy-exe
            KERNEL fl
            OPTIONS --friction-law=${law}
            LABELS proxy
        )
    endforeach()
    foreach(law 16 3tp 7)
        seissol_add_proxy_smoke_test(
            NAME "Proxy (clusters, friction law ${law}): all_dr"
            TARGET seissol-proxy-exe
            KERNEL all_dr
            OPTIONS --mode=cluster --friction-law=${law}
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
    # loops that do not go together, see seissol::proxy::checkConfig
    seissol_add_smoke_test(
        NAME "Proxy: a kernel of the kernel mode in the cluster mode"
        TARGET seissol-proxy-exe
        ARGS 10 1 ader --mode=cluster
        EXPECT_FAILURE
        LABELS proxy
    )
    seissol_add_smoke_test(
        NAME "Proxy: both kinds of plastic cells"
        TARGET seissol-proxy-exe
        ARGS 10 1 plast,plast_yield
        EXPECT_FAILURE
        LABELS proxy
    )
    seissol_add_smoke_test(
        NAME "Proxy: the plasticity of the cluster mode in the kernel mode"
        TARGET seissol-proxy-exe
        ARGS 10 1 all --plasticity=yield
        EXPECT_FAILURE
        LABELS proxy
    )
    seissol_add_smoke_test(
        NAME "Proxy: a friction law without parameters in the proxy"
        TARGET seissol-proxy-exe
        ARGS 10 1 fl --friction-law=2
        EXPECT_FAILURE
        LABELS proxy
    )
    seissol_add_smoke_test(
        NAME "Proxy: thermal pressurization without rate and state"
        TARGET seissol-proxy-exe
        ARGS 10 1 fl --friction-law=16tp
        EXPECT_FAILURE
        LABELS proxy
    )
endif()
