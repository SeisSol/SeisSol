# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Instruments the C and C++ code of the build with the sanitizers that SANITIZERS names: every
# translation unit of the targets defined after this file is included (SeisSol, its submodules, the
# generated kernels, the tests and the proxy), and every executable and shared library links the
# runtime. Device code (CUDA, HIP and SYCL kernels) is not instrumented.
#
# The sanitizers are configured at runtime, with ASAN_OPTIONS, UBSAN_OPTIONS, LSAN_OPTIONS and
# TSAN_OPTIONS; .ci/sanitizers holds the options and suppressions the CI runs with.

if (NOT SANITIZERS)
  return()
endif()

# as a list, separated by semicolons or commas
string(REPLACE "," ";" _sanitizers "${SANITIZERS}")

foreach(_sanitizer IN LISTS _sanitizers)
  if (NOT _sanitizer IN_LIST SANITIZERS_OPTIONS)
    message(FATAL_ERROR
      "SANITIZERS contains \"${_sanitizer}\", which is no sanitizer. Allowed: ${SANITIZERS_OPTIONS}")
  endif()
endforeach()

# the runtimes of these pairs cannot be linked into one executable
foreach(_pair address:thread address:memory leak:thread leak:memory thread:memory)
  string(REPLACE ":" ";" _pair "${_pair}")
  list(GET _pair 0 _first)
  list(GET _pair 1 _second)
  if ((_first IN_LIST _sanitizers) AND (_second IN_LIST _sanitizers))
    message(FATAL_ERROR "The sanitizers ${_first} and ${_second} cannot be combined.")
  endif()
endforeach()

if ("thread" IN_LIST _sanitizers AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # libgomp synchronizes its threads in a way the thread sanitizer does not see, so that it reports
  # races in every parallel region; the OpenMP runtime of LLVM tells the sanitizer about its
  # synchronization through Archer (libarcher)
  message(WARNING "The thread sanitizer reports races in every OpenMP region with libgomp. "
                  "Build with Clang and the OpenMP runtime of LLVM instead.")
endif()

if ("memory" IN_LIST _sanitizers)
  message(WARNING "The memory sanitizer reports reads of uninitialized memory wherever "
                  "uninstrumented code (the C++ standard library, MPI, HDF5, ...) wrote it. "
                  "All dependencies need to be built with SANITIZERS=memory as well.")
endif()

list(JOIN _sanitizers "," _sanitizer_list)
set(SEISSOL_SANITIZER_FLAGS -fsanitize=${_sanitizer_list} -fno-omit-frame-pointer)
if ("memory" IN_LIST _sanitizers)
  list(APPEND SEISSOL_SANITIZER_FLAGS -fsanitize-memory-track-origins)
endif()

# fail here rather than at the first link of the build if the compiler lacks a runtime
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
string(REPLACE "," "_" _sanitizer_check "SEISSOL_SANITIZERS_LINK_${_sanitizer_list}")
cmake_push_check_state(RESET)
list(JOIN SEISSOL_SANITIZER_FLAGS " " CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LINK_OPTIONS ${SEISSOL_SANITIZER_FLAGS})
check_cxx_source_compiles("int main() { return 0; }" ${_sanitizer_check})
cmake_pop_check_state()
if (NOT ${_sanitizer_check})
  message(FATAL_ERROR "${CMAKE_CXX_COMPILER} cannot build a program with ${SEISSOL_SANITIZER_FLAGS}; "
                      "the runtime of a sanitizer may be missing.")
endif()

message(STATUS "Sanitizers: ${_sanitizer_list}")

add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:${SEISSOL_SANITIZER_FLAGS}>")
add_link_options("$<$<LINK_LANGUAGE:C,CXX>:${SEISSOL_SANITIZER_FLAGS}>")
