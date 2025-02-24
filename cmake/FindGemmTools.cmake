# SPDX-FileCopyrightText: 2020-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#  GemmTools - BLAS-like Library Instantiation Software Framework
#  email: ravil.dorozhinskii@tum.de 
#
#  Input:
#  GEMM_TOOLS_LIST                - a list of targeted BLAS implementations
#                                   which support GEneral Matrix Matrix multiplications
#
#  Output:
#  GemmTools_FOUND                - system has found GemmTools
#  GemmTools_INCLUDE_DIRS         - include directories for GemmTools
#  GemmTools_LIBRARIES            - libraries for GemmTools
#  GemmTools_COMPILER_DEFINITIONS - compiler definitions for GemmTools 
#
#  Example usage:
#
#  set(GEMM_TOOLS_LIST "LIBXSMM,PSpaMM")
#
#  find_package(GEMM_TOOLS_LIST)
#  if(GEMM_TOOLS_FOUND)
#    if(GEMM_TOOLS_INCLUDE_DIRS)
#      target_include_directories(TARGET PUBLIC ${GEMM_TOOLS_INCLUDE_DIRS})
#    endif()
#    if(GEMM_TOOLS_LIBRARIES)
#      target_link_libraries(TARGET PUBLIC ${GEMM_TOOLS_LIBRARIES})
#    endif()
#    if(GEMM_TOOLS_COMPILER_DEFINITIONS)  
#      target_compile_definitions(TARGET PUBLIC ${GEMM_TOOLS_COMPILER_DEFINITIONS})
#    endif()
#  endif()

string(REPLACE "," ";" _GEMM_TOOLS_LIST ${GEMM_TOOLS_LIST})

set(GemmTools_INCLUDE_DIRS "")
set(GemmTools_LIBRARIES "")
set(GemmTools_COMPILER_DEFINITIONS "")

foreach(component ${_GEMM_TOOLS_LIST})
    if ("${component}" STREQUAL "LIBXSMM")
        find_package(Libxsmm_executable REQUIRED)

    elseif("${component}" STREQUAL "LIBXSMM_JIT")
        find_package(LIBXSMM 1.17 REQUIRED)
        find_package(BLAS REQUIRED)

        set(GemmTools_INCLUDE_DIRS ${GemmTools_INCLUDE_DIRS} ${LIBXSMM_INCLUDE_DIRS})
        set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} ${LIBXSMM_LIBRARIES} ${BLAS_LIBRARIES})

    elseif ("${component}" STREQUAL "PSpaMM")
        find_package(PSpaMM REQUIRED)

    elseif ("${component}" STREQUAL "MKL")
        # cf. https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2025-0/cmake-config-for-onemkl.html
        find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})

        set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} MKL::MKL)

    elseif ("${component}" STREQUAL "OpenBLAS")
        find_package(OpenBLAS REQUIRED)
        set(GemmTools_INCLUDE_DIRS ${GemmTools_INCLUDE_DIRS} ${OpenBLAS_INCLUDE_DIRS})
        set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} ${OpenBLAS_LIBRARIES} ${BLAS_LIBRARIES})

    elseif ("${component}" STREQUAL "BLIS")
        find_package(BLIS REQUIRED)
        set(GemmTools_INCLUDE_DIRS ${GemmTools_INCLUDE_DIRS} ${BLIS_INCLUDE_DIRS})
        set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} ${BLIS_LIBRARIES})

    elseif ("${component}" STREQUAL "Eigen")
        # already included by default!

    elseif ("${component}" MATCHES "^[ \t\r\n]*[Nn][Oo][Nn][Ee][ \t\r\n]*$")
        # (cf. https://cmake.org/cmake/help/latest/command/string.html#regex-specification)

        # no includes necessary

    elseif ("${component}" STREQUAL "GemmForge")
        execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import gemmforge; gemmforge.print_cmake_path()"
                        OUTPUT_VARIABLE GEMMFORGE_PATH)
        set(CMAKE_PREFIX_PATH "${GEMMFORGE_PATH}" ${CMAKE_PREFIX_PATH})
        find_package(GemmForge 0.0.207 REQUIRED)
        set(DEVICE_SRC ${DEVICE_SRC} ${GemmForge_SOURCES})
        set(DEVICE_INCLUDE_DIRS ${DEVICE_INCLUDE_DIRS} ${GemmForge_INCLUDE_DIRS})

    elseif ("${component}" STREQUAL "tinytc")
        find_package(tinytc 0.3.1 REQUIRED shared)
        find_package(tinytc_sycl 0.3.1 REQUIRED shared)
        set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} tinytc::tinytc tinytc::tinytc_sycl)

    elseif ("${component}" STREQUAL "TensorForge")
        execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import tensorforge; tensorforge.print_cmake_path()"
                        OUTPUT_VARIABLE TENSORFORGE_PATH)
        set(CMAKE_PREFIX_PATH "${TENSORFORGE_PATH}" ${CMAKE_PREFIX_PATH})
        find_package(TensorForge REQUIRED)
        set(DEVICE_SRC ${DEVICE_SRC} ${TensorForge_SOURCES})
        set(DEVICE_INCLUDE_DIRS ${DEVICE_INCLUDE_DIRS} ${TensorForge_INCLUDE_DIRS})

    else()
        message(FATAL_ERROR "Gemm Tools do not have a requested component, i.e. ${component}. \
                Please, refer to the documentation")
    endif()

endforeach()
