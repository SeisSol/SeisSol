# modified from https://gitlab.kitware.com/cmake/cmake/blob/master/Modules/FindCxxTest.cmake
# Original Copyright 2000-2020 Kitware, Inc. and Contributors
# SPDX-License-Identifier: BSD-3-Clause

macro(CXXTEST_ADD_TEST_MPI _cxxtest_testname _ranks _cxxtest_outfname)
    set(_cxxtest_real_outfname ${CMAKE_CURRENT_BINARY_DIR}/${_cxxtest_outfname})
    if (MPI)
        set(_cxxtest_command_prefix mpirun -n ${_ranks} --oversubscribe)
    else()
        set(_cxxtest_command_prefix)
    endif()

    add_custom_command(
            OUTPUT  ${_cxxtest_real_outfname}
            DEPENDS ${ARGN}
            COMMAND ${CXXTEST_TESTGEN_INTERPRETER}
            ${CXXTEST_TESTGEN_EXECUTABLE} ${CXXTEST_TESTGEN_ARGS} -o ${_cxxtest_real_outfname} ${ARGN}
    )

    set_source_files_properties(${_cxxtest_real_outfname} PROPERTIES GENERATED true)
    add_executable(${_cxxtest_testname} ${_cxxtest_real_outfname} ${ARGN})

    if(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
        add_test(
                NAME ${_cxxtest_testname}
                COMMAND {_cxxtest_command_prefix} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_cxxtest_testname}"
        )
    elseif(EXECUTABLE_OUTPUT_PATH)
        add_test(
                NAME ${_cxxtest_testname}
                COMMAND ${_cxxtest_command_prefix} "${EXECUTABLE_OUTPUT_PATH}/${_cxxtest_testname}")
    else()
        add_test(
                NAME ${_cxxtest_testname}
                COMMAND ${_cxxtest_command_prefix} "${CMAKE_CURRENT_BINARY_DIR}/${_cxxtest_testname}")
    endif()

endmacro()
