# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# ---------------------------------------------------------------------------
# Test data files
# ---------------------------------------------------------------------------

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Initializer/TimeStepping/mesh.h5
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/mesh.h5
    COPYONLY
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Initializer/TimeStepping/material.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/material.yaml
    COPYONLY
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/ResultWriter/receiver_correct.dat
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/receiver_correct.dat
    COPYONLY
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Reader/fsrm_source1.dat
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/fsrm_source1.dat
    COPYONLY
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Reader/fsrm_source2.dat
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/fsrm_source2.dat
    COPYONLY
)

if (NETCDF)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/Reader/source_loh.nrf
        ${CMAKE_CURRENT_BINARY_DIR}/Testing/source_loh.nrf
        COPYONLY
    )
endif()

configure_file("Setup.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/Setup.h")


# =============================================================================
# Shared test infrastructure
# =============================================================================

# Interface library that carries the common compile/link settings for all test
# targets.  Every per-module binary links against this instead of duplicating
# the same target_link_libraries / target_include_directories calls.
add_library(seissol-test-base INTERFACE)
target_link_libraries(seissol-test-base INTERFACE seissol-lib seissol-external)
target_include_directories(seissol-test-base INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(seissol-test-base INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(seissol-test-base INTERFACE YATETO_TESTING_NO_FLOP_COUNTER)

separate_arguments(TESTING_COMMAND_LIST NATIVE_COMMAND ${TESTING_COMMAND})

# =============================================================================
# Helper function: add a per-module test executable
# =============================================================================
#
# Usage:
#   seissol_add_test_module(<name> <source1> [<source2> ...])
#
# Creates an executable  seissol-test-<name>  from the given sources plus the
# shared TestMain.cpp.  The target is automatically registered with CTest via
# doctest_discover_tests().  Doctest test_suite labels are forwarded as CTest
# labels (ADD_LABELS) so that  ctest -L <suite>  works out of the box.
#

function(seissol_add_test_module NAME)
    set(TARGET_NAME "seissol-test-${NAME}")

    add_executable(${TARGET_NAME} ${ARGN} TestMain.cpp)
    target_link_libraries(${TARGET_NAME} PRIVATE seissol-test-base)

    separate_arguments(_CMD_LIST NATIVE_COMMAND ${TESTING_COMMAND})
    set_property(TARGET ${TARGET_NAME} PROPERTY CROSSCOMPILING_EMULATOR ${_CMD_LIST})

    if (NOT TESTING_BUILD_ONLY)
        doctest_discover_tests(${TARGET_NAME}
            ADD_LABELS 1
            TEST_PREFIX "${NAME}/"
        )
    endif()
endfunction()


# =============================================================================
# Helper function: add a per-module MPI test executable
# =============================================================================
#
# Usage:
#   seissol_add_mpi_test_module(<n> <ranks> <source1> [<source2> ...])
#
# Like seissol_add_test_module, but both discovery and execution run under
# mpirun -n <ranks>.  doctest_discover_tests tolerates the duplicate stdout
# from multiple ranks during --list-test-cases.
#
function(seissol_add_mpi_test_module NAME RANKS)
    set(TARGET_NAME "seissol-test-${NAME}")

    add_executable(${TARGET_NAME} ${ARGN} TestMain.cpp)
    target_link_libraries(${TARGET_NAME} PRIVATE seissol-test-base)

    set_property(TARGET ${TARGET_NAME}
        PROPERTY CROSSCOMPILING_EMULATOR
        mpirun -n ${RANKS} --oversubscribe --allow-run-as-root
    )

    if (NOT TESTING_BUILD_ONLY)
        doctest_discover_tests(${TARGET_NAME}
            ADD_LABELS 1
            TEST_PREFIX "${NAME}/"
        )
    endif()
endfunction()

# ---------------------------------------------------------------------------
# Per-module test targets
# ---------------------------------------------------------------------------

seissol_add_test_module(common          Common/TestCommon.cpp)
seissol_add_test_module(dynamicrupture  DynamicRupture/TestDynamicRupture.cpp)
seissol_add_test_module(equations       Equations/TestEquations.cpp)
seissol_add_test_module(geometry        Geometry/TestGeometry.cpp)
seissol_add_test_module(initializer     Initializer/TestInitializer.cpp)
seissol_add_test_module(io              IO/TestIO.cpp)
seissol_add_test_module(kernel          Kernel/TestKernel.cpp)
seissol_add_test_module(memory          Memory/TestMemory.cpp)
seissol_add_test_module(model           Model/TestModel.cpp)
seissol_add_test_module(modules         Modules/TestModule.cpp)
seissol_add_test_module(monitoring      Monitoring/TestMonitoring.cpp)
seissol_add_test_module(numerical       Numerical/TestNumerical.cpp)
seissol_add_test_module(parallel        Parallel/TestParallel.cpp)
seissol_add_test_module(physics         Physics/InstantaneousTimeMirror.cpp)
seissol_add_test_module(resultwriter    ResultWriter/TestResultWriter.cpp)
seissol_add_test_module(solver          Solver/TestSolver.cpp)
seissol_add_test_module(sourceterm      SourceTerm/TestSourceTerm.cpp)

if (NETCDF)
    seissol_add_test_module(reader      Reader/ReaderTest.cpp)
endif()

# =============================================================================
# Multi-rank MPI test targets (WS2)
# =============================================================================
#
# These tests are executed via  mpirun -n <N>  and exercise MPI collectives,
# parallel statistics, and cross-rank data redistribution.

seissol_add_mpi_test_module(mpi 4
    Parallel/TestParallelMPI.cpp
)

# =============================================================================
# Handle TESTING_GENERATED kernel sources
# =============================================================================
#
# The generated kernel test code must be compiled into every target that
# contains tests exercising codegen kernels.  GPU builds additionally need a
# separate device library.

if (TESTING_GENERATED)
    set(_GENERATED_KERNEL_TARGETS seissol-test-kernel)

    foreach(_TGT IN LISTS _GENERATED_KERNEL_TARGETS)
        if (WITH_GPU)
            make_device_lib(${_TGT}-device "${CMAKE_BINARY_DIR}/codegen/GeneratedCode/test-kernel.cpp")
            target_link_libraries(${_TGT}-device PRIVATE seissol-kernel-lib)
            target_link_libraries(${_TGT}-device PRIVATE seissol-external)
            target_include_directories(${_TGT}-device PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
            target_include_directories(${_TGT}-device PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
            target_link_libraries(${_TGT} PRIVATE ${_TGT}-device)
            target_compile_definitions(${_TGT}-device PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)
        else()
            target_sources(${_TGT} PRIVATE "${CMAKE_BINARY_DIR}/codegen/GeneratedCode/test-kernel.cpp")
        endif()
    endforeach()
endif()

# Get all propreties that cmake supports
if(NOT CMAKE_PROPERTY_LIST)
    execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

    # Convert command output into a CMake list
    string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)
endif()

if (TESTING_PYTHON)
    add_subdirectory(python/codegen)
endif()
