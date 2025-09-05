# SPDX-FileCopyrightText: 2021 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

set(DEVICE_SRC ${DEVICE_SRC}
          ${CMAKE_BINARY_DIR}/codegen/GeneratedCode/gpulike_subroutine.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/sycl/PlasticityAux.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/LinearCK/DeviceAux/sycl/KernelsAux.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverSycl.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/PointSourceClusterSycl.cpp)

if (("${DEVICE_BACKEND}" STREQUAL "hipsycl") OR ("${DEVICE_BACKEND}" STREQUAL "acpp"))

    find_package(Boost REQUIRED COMPONENTS context fiber)
    if (DEVICE_ARCH MATCHES "sm_*")
        if (CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
            set(SYCL_USE_NVHPC_DEFAULT ON)
        else()
            set(SYCL_USE_NVHPC_DEFAULT OFF)
        endif()
        option(SYCL_USE_NVHPC "For Nvidia GPUs, use nvhpc instead of CUDA/nvcc." ${SYCL_USE_NVHPC_DEFAULT})
        if (SYCL_USE_NVHPC)
            # we assume that hipsycl was compiled with nvhpc compiler collection
            string(REPLACE sm_ cc NVCPP_ARCH ${DEVICE_ARCH})
            set(HIPSYCL_TARGETS "cuda-nvcxx:${NVCPP_ARCH}" CACHE STRING "" FORCE)
            set(ACPP_TARGETS "cuda-nvcxx:${NVCPP_ARCH}" CACHE STRING "" FORCE)
        else()
            set(HIPSYCL_TARGETS "cuda:${DEVICE_ARCH}" CACHE STRING "" FORCE)
            set(ACPP_TARGETS "cuda:${DEVICE_ARCH}" CACHE STRING "" FORCE)
            target_compile_options(device PRIVATE -Wno-unknown-cuda-version)
        endif()
    elseif(DEVICE_ARCH MATCHES "gfx*")
        set(HIPSYCL_TARGETS "hip:${DEVICE_ARCH}" CACHE STRING "" FORCE)
        set(ACPP_TARGETS "hip:${DEVICE_ARCH}" CACHE STRING "" FORCE)
    else()
        set(HIPSYCL_TARGETS "${DEVICE_BACKEND}:${DEVICE_ARCH}" CACHE STRING "" FORCE)
        set(ACPP_TARGETS "${DEVICE_BACKEND}:${DEVICE_ARCH}" CACHE STRING "" FORCE)
    endif()

    find_package(OpenMP REQUIRED)
    find_package(AdaptiveCpp CONFIG REQUIRED)

    function(make_device_lib NAME FILES)
        add_library(${NAME} SHARED ${FILES})
        target_include_directories(${NAME} PUBLIC ${SEISSOL_DEVICE_INCLUDE})

        target_link_libraries(${NAME} PUBLIC ${Boost_LIBRARIES})

        target_compile_options(${NAME} PRIVATE ${EXTRA_CXX_FLAGS} "-fPIC" ${OpenMP_CXX_FLAGS})

        add_sycl_to_target(TARGET ${NAME} SOURCES ${FILES})
    endfunction()

elseif("${DEVICE_BACKEND}" STREQUAL "oneapi")
    find_package(DpcppFlags REQUIRED)

    function(make_device_lib NAME FILES)
        add_library(${NAME} SHARED ${FILES})
        target_include_directories(${NAME} PUBLIC ${SEISSOL_DEVICE_INCLUDE})
        target_compile_options(${NAME} PRIVATE ${EXTRA_CXX_FLAGS} "-O3")
        target_compile_definitions(${NAME} PRIVATE __DPCPP_COMPILER)
        target_link_libraries(${NAME} PRIVATE dpcpp::device_flags ${GemmTools_LIBRARIES})
    endfunction()
    target_link_libraries(seissol-common-properties INTERFACE dpcpp::interface)
endif()

make_device_lib(seissol-device-lib "${DEVICE_SRC}")
