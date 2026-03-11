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

# ---------------------------------------------------------------------------
# Shared test infrastructure (INTERFACE – headers only, no object code)
# ---------------------------------------------------------------------------

add_library(seissol-test-base INTERFACE)
target_link_libraries(seissol-test-base INTERFACE seissol-lib seissol-external)
target_include_directories(seissol-test-base INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
)

separate_arguments(TESTING_COMMAND_LIST NATIVE_COMMAND ${TESTING_COMMAND})

# ---------------------------------------------------------------------------
# Macro: create a per-module test executable
# ---------------------------------------------------------------------------
# Usage:   seissol_add_test_module(numerical  Numerical/TestNumerical.cpp)
#
# Creates target  seissol-test-<name>  that is independently buildable,
# runnable, and discoverable by ctest.
# ---------------------------------------------------------------------------

macro(seissol_add_test_module _name)
    set(_target "seissol-test-${_name}")
    add_executable(${_target} ${ARGN} TestMain.cpp)
    target_link_libraries(${_target} PRIVATE seissol-test-base)
    target_compile_definitions(${_target} PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)
    set_property(TARGET ${_target} PROPERTY CROSSCOMPILING_EMULATOR ${TESTING_COMMAND_LIST})
    if (NOT TESTING_BUILD_ONLY)
        doctest_discover_tests(${_target})
    endif()
endmacro()

# ---------------------------------------------------------------------------
# Per-module test targets
# ---------------------------------------------------------------------------

seissol_add_test_module(common          Common/TestCommon.cpp)
seissol_add_test_module(dynamicrupture  DynamicRupture/TestDynamicRupture.cpp)
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
seissol_add_test_module(resultwriter    ResultWriter/TestResultWriter.cpp)
seissol_add_test_module(solver          Solver/TimeStepping/TestSolverTimeStepping.cpp)
seissol_add_test_module(sourceterm      SourceTerm/TestSourceTerm.cpp)

if (NETCDF)
    seissol_add_test_module(reader      Reader/ReaderTest.cpp)
endif()

# ---------------------------------------------------------------------------
# Monolithic binary  (almost obsolete; keep for yateto tests for now)
# ---------------------------------------------------------------------------
# Kept so that the Yateto tests continue to work.
# Once those are also added to the CI, this
# section can be removed.
# ---------------------------------------------------------------------------

set(seissol_test_sources
)

add_executable(seissol-serial-test
        ${seissol_test_sources}
        TestMain.cpp)
target_link_libraries(seissol-serial-test PRIVATE seissol-lib)
target_link_libraries(seissol-serial-test PUBLIC seissol-external)
target_include_directories(seissol-serial-test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(seissol-serial-test PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if (TESTING_GENERATED)
    if (WITH_GPU)

        # for the GPU, create an extra library with device linking flags etc.
        make_device_lib(seissol-serial-test-device "${CMAKE_BINARY_DIR}/codegen/GeneratedCode/test-kernel.cpp")
        target_link_libraries(seissol-serial-test-device PRIVATE seissol-kernel-lib)
        target_link_libraries(seissol-serial-test-device PRIVATE seissol-external)
        target_include_directories(seissol-serial-test-device PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
        target_include_directories(seissol-serial-test-device PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
        target_link_libraries(seissol-serial-test PRIVATE seissol-serial-test-device)
        target_compile_definitions(seissol-serial-test-device PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)

    else()

        # else (i.e. CPU only), just add it to the main test binary
        target_sources(seissol-serial-test PRIVATE "${CMAKE_BINARY_DIR}/codegen/GeneratedCode/test-kernel.cpp")
    endif()
endif()

set_property(TARGET seissol-serial-test PROPERTY CROSSCOMPILING_EMULATOR ${TESTING_COMMAND_LIST})
if (NOT TESTING_BUILD_ONLY)
    doctest_discover_tests(seissol-serial-test)
endif()

# Avoid duplicate definition of FLOP counters
target_compile_definitions(seissol-serial-test PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)

# Get all propreties that cmake supports
if(NOT CMAKE_PROPERTY_LIST)
    execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

    # Convert command output into a CMake list
    string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)
endif()
