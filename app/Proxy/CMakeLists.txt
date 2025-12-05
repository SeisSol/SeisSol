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
    add_test(NAME "Proxy: all" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 all)
    add_test(NAME "Proxy: all_dr" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 all_dr)
    add_test(NAME "Proxy: ader" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 ader)
    add_test(NAME "Proxy: localwoader" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 localwoader)
    add_test(NAME "Proxy: local" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 local)
    add_test(NAME "Proxy: neigh" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 neigh)
    add_test(NAME "Proxy: neigh_dr" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 neigh_dr)
    add_test(NAME "Proxy: godunov_dr" COMMAND $<TARGET_FILE:seissol-proxy-exe> 10 1 godunov_dr)
endif()
