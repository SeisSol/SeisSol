# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

set(seissol_test_sources
        Model/TestModel.cpp
        Initializer/TestInitializer.cpp
        Numerical/TestNumerical.cpp
        Geometry/TestGeometry.cpp
        Kernel/TestKernel.cpp
        Memory/TestMemory.cpp
        SourceTerm/TestSourceTerm.cpp
        Parallel/TestParallel.cpp
        ResultWriter/TestResultWriter.cpp
        Solver/TimeStepping/TestSolverTimeStepping.cpp
        DynamicRupture/TestDynamicRupture.cpp
        Common/TestCommon.cpp
)

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
    set(seissol_test_sources ${seissol_test_sources} ${CMAKE_CURRENT_SOURCE_DIR}/Reader/ReaderTest.cpp)
endif()

add_executable(seissol-serial-test
        ${seissol_test_sources}
        TestMain.cpp)
target_link_libraries(seissol-serial-test PRIVATE seissol-lib)
target_link_libraries(seissol-serial-test PUBLIC seissol-external)
target_include_directories(seissol-serial-test PUBLIC .)

if (TESTING_GENERATED)
    if (WITH_GPU)

        # for the GPU, create an extra library with device linking flags etc.
        make_device_lib(seissol-serial-test-device "${CMAKE_BINARY_DIR}/codegen/GeneratedCode/test-kernel.cpp")
        target_link_libraries(seissol-serial-test-device PRIVATE seissol-kernel-lib)
        target_link_libraries(seissol-serial-test-device PRIVATE seissol-external)
        target_include_directories(seissol-serial-test-device PRIVATE .)
        target_link_libraries(seissol-serial-test PRIVATE seissol-serial-test-device)
        target_compile_definitions(seissol-serial-test-device PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)

    else()

        # else (i.e. CPU only), just add it to the main test binary
        target_sources(seissol-serial-test PRIVATE "${CMAKE_BINARY_DIR}/codegen/GeneratedCode/test-kernel.cpp")
    endif()
endif()

separate_arguments(TESTING_COMMAND_LIST NATIVE_COMMAND ${TESTING_COMMAND})
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
