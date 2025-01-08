# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

add_subdirectory(Proxy)

if (PROXY_PYBINDING)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        find_package(pybind11 2.6.2 REQUIRED)
    else()
        find_package(pybind11 2.3.0 REQUIRED)
    endif()

    pybind11_add_module(seissol_proxy_bindings Bindings.cpp)
    target_link_libraries(seissol_proxy_bindings PUBLIC seissol-proxy-lib)
    file(COPY proxy-runners DESTINATION ${CMAKE_BINARY_DIR})
endif()

if (WITH_GPU)
    target_link_libraries(seissol-proxy-lib PUBLIC general-sycl-offloading)
endif()

# C/C++ proxy interface
add_executable(seissol-proxy Main.cpp)
target_link_libraries(seissol-proxy PUBLIC seissol-proxy-lib)
set_target_properties(seissol-proxy PROPERTIES OUTPUT_NAME "SeisSol_proxy_${EXE_NAME_PREFIX}")
install(TARGETS seissol-proxy RUNTIME)
