#  OpenMP_Offloading - finds flags for OpenMP Offloading
#
#  requires 1. find_package(OpenMP REQUIRED) call before
#           2. DEVICE_BACKEND must be set
#           3. DEVICE_ARCH must be set
#
#  OpenMP_Offloading_FOUND               - compiler is able to perform OpenMP offloading
#  OpenMP_Offloading_FLAGS               - compiler flags for OpenMP Offloading
#  OpenMP_Offloading_LINK_LIBRARIES      - link libraries for OpenMP Offloading
#  OpenMP::OpenMP_CXX_with_offloading    - imported target
#
#  NOTE: it is not an official cmake search file
#  authors: Ravil Dorozhinskii
#  email: ravil.dorozhinskii@tum.de
#

include(FindPackageHandleStandardArgs)
include(CheckCXXSourceRuns)
include(CMakePushCheckState)

cmake_push_check_state()

set(_is_offloading_test_ready True)
if (NOT OpenMP_FOUND)
    message(WARNING "OpenMP is not found. Make sure that you called find_package(OpenMP REQUIRED) first")
    set(_is_offloading_test_ready False)
endif()

if (NOT DEFINED DEVICE_BACKEND)
    message(WARNING "DEVICE_BACKEND is not defined")
    set(_is_offloading_test_ready False)
endif()

if (NOT DEFINED DEVICE_ARCH)
    message(WARNING "DEVICE_ARCH is not defined")
    set(_is_offloading_test_ready False)
endif()

if (_is_offloading_test_ready)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        if (${DEVICE_BACKEND} STREQUAL "cuda")
            set(_Offloading_CXX_FLAGS -fopenmp-targets=nvptx64-nvidia-cuda
                                      -Xopenmp-target
                                      -march=${DEVICE_ARCH})
            set(_Offloading_LINK_LIBS -lomptarget)
        endif()
    elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
        if (${DEVICE_BACKEND} STREQUAL "cuda")
            set(_Offloading_CXX_FLAGS --target=nvptx-none)
            set(_Offloading_LINK_LIBS -lomptarget)
        endif()
    endif()

    set(CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS} ${_Offloading_CXX_FLAGS})
    string (REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")

    set(CMAKE_REQUIRED_INCLUDES ${OpenMP_CXX_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_CXX_LIBRARIES})

    check_cxx_source_runs("
    #include <omp.h>
    #define YES 0
    #define NO -1
    int main(int argc, char* argv[]) {
      int canOffload = NO;
      #pragma omp target map(tofrom: canOffload)
      {
        if (!omp_is_initial_device()) {
          canOffload = YES;
        }
      }
      return canOffload;
    }
    " _OPENMP_OFFLOAD_RESULT)
else()
    message(WARNING "OpenMP Offloading test skipped")
endif()

if (_OPENMP_OFFLOAD_RESULT)
    set(OpenMP_Offloading_FLAGS ${_Offloading_CXX_FLAGS})
    set(OpenMP_Offloading_LINK_LIBRARIES ${_Offloading_LINK_LIBS})

    add_library(OpenMP::OpenMP_CXX_with_offloading INTERFACE IMPORTED)
    set_property(TARGET OpenMP::OpenMP_CXX_with_offloading APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
            $<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS} ${_Offloading_CXX_FLAGS}>)

    target_link_libraries(OpenMP::OpenMP_CXX_with_offloading INTERFACE
            ${OpenMP_CXX_LIBRARIES}
            ${OpenMP_CXX_FLAGS} ${OpenMP_Offloading_FLAGS}
            ${OpenMP_Offloading_LINK_LIBRARIES})
endif()

find_package_handle_standard_args(OpenMP_Offloading
        OpenMP_Offloading_FLAGS
        OpenMP_Offloading_LINK_LIBRARIES
        _OPENMP_OFFLOAD_RESULT)
mark_as_advanced(_is_offloading_test_ready _Offloading_CXX_FLAGS _Offloading_LINK_LIBS)
cmake_pop_check_state()
