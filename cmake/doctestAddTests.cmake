# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

# Imported copyright info from there:
# SPDX-FileCopyrightText: 2000-2025 Kitware, Inc. and Contributors

set(prefix "${TEST_PREFIX}")
set(suffix "${TEST_SUFFIX}")
set(spec ${TEST_SPEC})
set(extra_args ${TEST_EXTRA_ARGS})
set(properties ${TEST_PROPERTIES})
set(script)
set(suite)
set(tests)

function(add_command NAME)
  set(_args "")
  foreach(_arg ${ARGN})
    if(_arg MATCHES "[^-./:a-zA-Z0-9_]")
      set(_args "${_args} [==[${_arg}]==]") # form a bracket_argument
    else()
      set(_args "${_args} ${_arg}")
    endif()
  endforeach()
  set(script "${script}${NAME}(${_args})\n" PARENT_SCOPE)
endfunction()

# Run test executable to get list of available tests
if(NOT EXISTS "${TEST_EXECUTABLE}")
  message(FATAL_ERROR
    "Specified test executable '${TEST_EXECUTABLE}' does not exist"
  )
endif()

if("${spec}" MATCHES .)
  set(spec "--test-case=${spec}")
endif()

# In a build with sanitizers (cmake/sanitizers.cmake), listing the tests runs MPI_Init and
# MPI_Finalize, after which Open MPI leaves memory allocated, and in which it takes its locks in an
# order that the thread sanitizer reports without the options the tests run with
# (.ci/sanitizers/environment.sh); what the sanitizers find is for the runs of the tests to
# report, not for the listing of them.
foreach(_variable_and_option ASAN_OPTIONS:detect_leaks=0 LSAN_OPTIONS:detect_leaks=0
                             TSAN_OPTIONS:report_bugs=0)
  string(REPLACE ":" ";" _variable_and_option "${_variable_and_option}")
  list(GET _variable_and_option 0 _variable)
  list(GET _variable_and_option 1 _option)
  if(DEFINED ENV{${_variable}})
    set(ENV{${_variable}} "$ENV{${_variable}}:${_option}")
  else()
    set(ENV{${_variable}} "${_option}")
  endif()
endforeach()
# With the thread sanitizer, PMIx (of Open MPI) cannot attach its shared memory at the address it
# wants, which the sanitizer keeps for itself, and MPI_Init crashes; PMIx stores its data without
# shared memory then (as .ci/sanitizers/environment.sh has the runs of the tests do).
if(NOT DEFINED ENV{PMIX_MCA_gds})
  set(ENV{PMIX_MCA_gds} "hash")
endif()

execute_process(
  COMMAND ${TEST_EXECUTOR} "${TEST_EXECUTABLE}" ${spec} --list-test-cases -nc
  OUTPUT_VARIABLE output
  RESULT_VARIABLE result
)
if(NOT ${result} EQUAL 0)
  message(FATAL_ERROR
    "Error running test executable '${TEST_EXECUTABLE}':\n"
    "  Result: ${result}\n"
    "  Output: ${output}\n"
  )
endif()

string(REPLACE "\n" ";" output "${output}")

# Parse output
foreach(line ${output})
  if("${line}" STREQUAL "===============================================================================" OR "${line}" MATCHES [==[^\[doctest\] ]==])
    continue()
  endif()
  set(test ${line})
  # use escape commas to handle properly test cases with commas inside the name
  string(REPLACE "," "\\," test_name ${test})
  # ...and add to script
  add_command(add_test
    "${prefix}${test}${suffix}"
    ${TEST_EXECUTOR}
    "${TEST_EXECUTABLE}"
    "--test-case=${test_name}"
    ${extra_args}
  )
  add_command(set_tests_properties
    "${prefix}${test}${suffix}"
    PROPERTIES
    WORKING_DIRECTORY "${TEST_WORKING_DIR}"
    ${properties}
  )
  list(APPEND tests "${prefix}${test}${suffix}")
endforeach()

# Create a list of all discovered tests, which users may use to e.g. set
# properties on the tests
add_command(set ${TEST_LIST} ${tests})

# Write CTest script
file(WRITE "${CTEST_FILE}" "${script}")
