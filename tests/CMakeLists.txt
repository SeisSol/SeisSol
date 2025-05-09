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
        SourceTerm/TestSourceTerm.cpp
        Parallel/TestParallel.cpp
        ResultWriter/TestResultWriter.cpp
        Solver/time_stepping/TestSolverTimeStepping.cpp
        DynamicRupture/TestDynamicRupture.cpp
        Common/TestCommon.cpp
        )


if (TESTING_GENERATED)
    # FIXME: fix by removing the `using namespace` in the test files
    # and replace by ${CMAKE_CURRENT_BINARY_DIR}/generated-code/generated_code/test-kernel.cpp
    file(GLOB_RECURSE KERNEL_TESTFILES ${CMAKE_BINARY_DIR}/generated-code/generated_code/**/test-kernel.cpp)
    set(seissol_test_sources
        ${seissol_test_sources}
        ${KERNEL_TESTFILES}
        )
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Initializer/time_stepping/mesh.h5
    ${CMAKE_CURRENT_BINARY_DIR}/Testing/mesh.h5
    COPYONLY
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Initializer/time_stepping/material.yaml
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
target_include_directories(seissol-serial-test PUBLIC ${CMAKE_SOURCE_DIR}/external/)

separate_arguments(TESTING_COMMAND_LIST NATIVE_COMMAND ${TESTING_COMMAND})
set_property(TARGET seissol-serial-test PROPERTY CROSSCOMPILING_EMULATOR ${TESTING_COMMAND_LIST})
doctest_discover_tests(seissol-serial-test)

# Avoid duplicate definition of FLOP counters
target_compile_definitions(seissol-serial-test PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)
