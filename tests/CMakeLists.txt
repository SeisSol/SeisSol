# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

set(SEISSOL_TEST_SOURCES
    Model/TestModel.cpp
    Initializer/TestInitializer.cpp
    Numerical/TestNumerical.cpp
    Geometry/TestGeometry.cpp
    Kernel/TestKernel.cpp
    SourceTerm/TestSourceTerm.cpp
    Parallel/TestParallel.cpp
    Pipeline/TestPipeline.cpp
    ResultWriter/TestResultWriter.cpp
    Solver/TimeStepping/TestSolverTimeStepping.cpp
    DynamicRupture/TestDynamicRupture.cpp
    Common/TestCommon.cpp
    )


if (TESTING_GENERATED)
    set(SEISSOL_TEST_SOURCES
        ${SEISSOL_TEST_SOURCES}
        ${CMAKE_BINARY_DIR}/generated-code/generated-code/test-kernel.cpp
        )
endif()

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
    set(SEISSOL_TEST_SOURCES ${SEISSOL_TEST_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/Reader/ReaderTest.cpp)
endif()

add_executable(SeisSol-serial-test
        ${SEISSOL_TEST_SOURCES}
        TestMain.cpp)
target_link_libraries(SeisSol-serial-test PRIVATE SeisSol-lib)
if(WITH_GPU)
    target_link_libraries(SeisSol-serial-test PRIVATE general-sycl-offloading)
endif()
target_include_directories(SeisSol-serial-test PUBLIC ../external/)
target_include_directories(SeisSol-serial-test PRIVATE .)
doctest_discover_tests(SeisSol-serial-test)

# Avoid duplicate definition of FLOP counters
target_compile_definitions(SeisSol-serial-test PRIVATE YATETO_TESTING_NO_FLOP_COUNTER)
