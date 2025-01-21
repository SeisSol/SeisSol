# SPDX-FileCopyrightText: 2020-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#  FILESYSTEM - filesystem component of the standard library
#
#  FILESYSTEM_FOUND       - system has std::filesystem
#  FILESYSTEM_LIBRARIES   - libraries for std::filesystem
#  std::filesystem        - imported target
#
#  NOTE: it is not an official cmake search file
#  authors: Carsten Uphoff, Ravil Dorozhinskii
#  email: uphoff@in.tum.de, ravil.dorozhinskii@tum.de 
#


include(FindPackageHandleStandardArgs)
include(CheckCXXSourceRuns)
include(CMakePushCheckState)

cmake_push_check_state()

set(_PARENT_CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD})
set(CMAKE_CXX_STANDARD 17)


set(_FILESYSTEM_TEST_RPOGRAM "
#include <iostream>

#ifdef EXPERIMENTAL_FS
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

int main(int argc, char* argv[]) {
    fs::path p{argv[0]};
    std::cout << fs::canonical(p) << std::endl;
    return 0;
}")

function(test_program filesystem_libraries use_experimental_fs is_ok)
    cmake_push_check_state(RESET)
    set(${is_ok} TRUE PARENT_SCOPE)
    if (${use_experimental_fs})
        set(CMAKE_REQUIRED_DEFINITIONS "-DEXPERIMENTAL_FS")
        set(TEST_SUFFIX "EXPERIMENTAL")
    else()
        set(TEST_SUFFIX "NATIVE")
    endif()

    check_cxx_source_runs("${_FILESYSTEM_TEST_RPOGRAM}" _FILESYSTEM_${TEST_SUFFIX})

    if (_FILESYSTEM_${TEST_SUFFIX})
        # a compiler has a native support for std::filesystem
        set(${filesystem_libraries} "" PARENT_SCOPE)
    else()
        set(CMAKE_REQUIRED_LIBRARIES "stdc++fs")
        check_cxx_source_runs("${_FILESYSTEM_TEST_RPOGRAM}" _FILESYSTEM_STDCPPFS_${TEST_SUFFIX})

        if(_FILESYSTEM_STDCPPFS_${TEST_SUFFIX})
            set(${filesystem_libraries} "stdc++fs" PARENT_SCOPE)
        else()
            set(CMAKE_REQUIRED_LIBRARIES "c++fs")
            check_cxx_source_runs("${_FILESYSTEM_TEST_RPOGRAM}" _FILESYSTEM_CPPFS_${TEST_SUFFIX})

            if(_FILESYSTEM_CPPFS_${TEST_SUFFIX})
                set(${filesystem_libraries} "c++fs" PARENT_SCOPE)
            else()
                set(${is_ok} FALSE PARENT_SCOPE)
            endif()
        endif()
    endif()
    cmake_pop_check_state()
endfunction()

set(FILESYSTEM_LIBRARIES "")
set(IS_OK FALSE)
test_program(FILESYSTEM_LIBRARIES FALSE IS_OK)

IF (IS_OK)
    set(IS_STD_FILESYSTEM_SUPPORT TRUE)
else()
    set(FILESYSTEM_LIBRARIES "")
    set(USE_EXPERIMENTAL_FS TRUE)
    test_program(FILESYSTEM_LIBRARIES TRUE IS_OK)

    if (IS_OK)
        set(IS_STD_FILESYSTEM_SUPPORT FALSE)
    else()
        # there is no way to compile a sample program with std::filesystem
        message(FATAL_ERROR "C++ compiler does not support std[::experimental]::filesystem")
    endif()
endif()


if (DEFINED FILESYSTEM_LIBRARIES)
    add_library(std::filesystem INTERFACE IMPORTED)
    set_property(TARGET std::filesystem APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_17)

    if (FILESYSTEM_LIBRARIES)
        set_property(TARGET std::filesystem APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${FILESYSTEM_LIBRARIES})
    endif()

    if (NOT IS_STD_FILESYSTEM_SUPPORT)
        set_property(TARGET std::filesystem APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS EXPERIMENTAL_FS)
    endif()
endif()

find_package_handle_standard_args(FILESYSTEM FILESYSTEM_LIBRARIES)

mark_as_advanced(FILESYSTEM_LIBRARIES _PARENT_CMAKE_CXX_STANDARD)
set(CMAKE_CXX_STANDARD ${_PARENT_CMAKE_CXX_STANDARD})
cmake_pop_check_state()
