# Copyright (C) 2023 Intel Corporation
# SPDX-License-Identifier: BSD-3-Clause

find_package(MPI REQUIRED)

set(SOURCES
    Args.cpp
    JsonLexer.cpp
    JsonToken.cpp
    LtsStructure.cpp
    LtsStructureParser.cpp
    MemoryType.cpp
    Simulator.cpp
    main.cpp
)

set(SEISSOL_SOURCES
    seissol/Parallel/MPI.cpp
    seissol/Solver/time_stepping/AbstractTimeCluster.cpp
    seissol/Solver/time_stepping/AbstractGhostTimeCluster.cpp
    seissol/Solver/time_stepping/ActorState.cpp
    seissol/Solver/time_stepping/DirectGhostTimeCluster.cpp
)

add_executable(SeisSol_comm_proxy ${SOURCES} ${SEISSOL_SOURCES})
target_compile_features(SeisSol_comm_proxy PRIVATE cxx_std_17)
target_link_libraries(SeisSol_comm_proxy PRIVATE MPI::MPI_CXX)
target_include_directories(SeisSol_comm_proxy PRIVATE . seissol)
target_compile_definitions(SeisSol_comm_proxy PRIVATE USE_MPI REAL_SIZE=4 LOG_LEVEL=2)
if(SYCL)
    find_package(DpcppFlags REQUIRED)
    target_link_libraries(SeisSol_comm_proxy PRIVATE dpcpp::device_flags)
    target_compile_definitions(SeisSol_comm_proxy PRIVATE HAVE_SYCL)
endif()
