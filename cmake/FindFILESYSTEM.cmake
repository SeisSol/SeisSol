#  FILESYSTEM - filesystem component of the standard library
#
#  FILESYSTEM_FOUND       - system has std::filesystem
#  FILESYSTEM_LIBRARIES   - libraries for std::filesystem
#
#  NOTE: it is not an official cmake search file
#  authors: Carsten Uphoff, Ravil Dorozhinskii
#  email: uphoff@in.tum.de, ravil.dorozhinskii@tum.de 
#


include(FindPackageHandleStandardArgs)
include(CheckCXXSourceRuns)

set(_FILESYSTEM_TEST_RPOGRAM "
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    fs::path p{argv[0]};
    std::cout << fs::canonical(p) << std::endl;
    return 0;
}")

set(CMAKE_REQUIRED_DEFINITIONS -std=c++17)
check_cxx_source_runs("${_FILESYSTEM_TEST_RPOGRAM}" _FILESYSTEM_NATIVE)

if (NOT _FILESYSTEM_NATIVE)
    set(CMAKE_REQUIRED_LIBRARIES "stdc++fs")
    check_cxx_source_runs("${_FILESYSTEM_TEST_RPOGRAM}" _FILESYSTEM_STDCPPFS)

    if(_FILESYSTEM_STDCPPFS)
        set(FILESYSTEM_LIBRARIES "stdc++fs")
    else()
        set(CMAKE_REQUIRED_LIBRARIES "c++fs")
        check_cxx_source_runs("${_FILESYSTEM_TEST_RPOGRAM}" _FILESYSTEM_CPPFS_RESULT)

        if(_FILESYSTEM_CPPFS_RESULT)
            set(FILESYSTEM_LIBRARIES "c++fs")
        else()
            # there is no way to compile a sample program with std::filesystem
            message(FATAL_ERROR "C++ compiler does not support std::filesystem")
        endif()
    endif()
else()
    # a compiler has a native support for std::filesystem
    set(FILESYSTEM_LIBRARIES "")
endif()


find_package_handle_standard_args(FILESYSTEM FILESYSTEM_LIBRARIES)

mark_as_advanced(FILESYSTEM_LIBRARIES)
