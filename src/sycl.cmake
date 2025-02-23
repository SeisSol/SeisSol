# SPDX-FileCopyrightText: 2021-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

set(DEVICE_SRC ${DEVICE_SRC}
          ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/sycl/PlasticityAux.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DeviceAux/sycl/KernelsAux.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverSycl.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/PointSourceClusterSycl.cpp)

add_library(seissol-device-lib SHARED ${DEVICE_SRC})

target_include_directories(seissol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
target_link_libraries(seissol-device-lib PRIVATE seissol-common-properties)

if ("${DEVICE_BACKEND}" STREQUAL "hipsycl")
  add_library(seissol-device-lib SHARED ${DEVICE_SRC})

  find_package(Boost REQUIRED COMPONENTS context fiber)
  if (DEVICE_ARCH MATCHES "sm_*")
    find_package(CUDA REQUIRED)
    set(HIPSYCL_TARGETS "cuda:${DEVICE_ARCH}")
    target_include_directories(seissol-device-lib PRIVATE ${CUDA_TOOLKIT_ROOT_DIR})
    target_link_libraries(seissol-device-lib PRIVATE -lcuda ${CUDA_LIBRARIES})
  elseif(DEVICE_ARCH MATCHES "gfx*")
    set(HIP_COMPILER hcc)
    find_package(HIP REQUIRED)
    set(HIPSYCL_TARGETS "hip:${DEVICE_ARCH}")
  else()
    set(HIPSYCL_TARGETS "${DEVICE_BACKEND}:${DEVICE_ARCH}")
  endif()

  target_compile_definitions(seissol-device-lib PRIVATE DEVICE_HIPSYCL_LANG REAL_SIZE=${REAL_SIZE_IN_BYTES})
  target_link_libraries(seissol-device-lib PUBLIC ${Boost_LIBRARIES})

  find_package(OpenMP REQUIRED)
  target_compile_options(seissol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-fPIC" ${OpenMP_CXX_FLAGS})

  find_package(hipSYCL CONFIG REQUIRED)
  add_sycl_to_target(TARGET seissol-device-lib SOURCES ${DEVICE_SRC})

elseif("${DEVICE_BACKEND}" STREQUAL "oneapi")
  find_package(DpcppFlags REQUIRED)
  target_compile_options(seissol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-O3")
  target_compile_definitions(seissol-device-lib PRIVATE DEVICE_ONEAPI_LANG REAL_SIZE=${REAL_SIZE_IN_BYTES} __DPCPP_COMPILER)
  target_link_libraries(seissol-device-lib PRIVATE dpcpp::device_flags ${GemmTools_LIBRARIES})
  target_link_libraries(seissol-common-properties INTERFACE dpcpp::interface)
endif()

target_compile_definitions(seissol-device-lib PRIVATE ${HARDWARE_DEFINITIONS}
        CONVERGENCE_ORDER=${ORDER}
        NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
        NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
        ${DR_QUAD_RULE})

if (DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
target_compile_definitions(seissol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
endif()
