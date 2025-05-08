# SPDX-FileCopyrightText: 2021 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

set(DEVICE_SRC ${DEVICE_SRC}
          ${CMAKE_BINARY_DIR}/generated-code/generated_code/gpulike_subroutine.cpp
          Kernels/DeviceAux/sycl/PlasticityAux.cpp
          Kernels/LinearCK/DeviceAux/sycl/KernelsAux.cpp
          DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverSycl.cpp
          Kernels/PointSourceClusterSycl.cpp)

add_library(seissol-device-lib SHARED ${DEVICE_SRC})

target_include_directories(seissol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
target_link_libraries(seissol-device-lib PRIVATE seissol-common-properties)

if (("${DEVICE_BACKEND}" STREQUAL "hipsycl") OR ("${DEVICE_BACKEND}" STREQUAL "acpp"))
  add_library(seissol-device-lib SHARED ${DEVICE_SRC})

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

  target_link_libraries(seissol-device-lib PUBLIC ${Boost_LIBRARIES})

  find_package(OpenMP REQUIRED)
  target_compile_options(seissol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-fPIC" ${OpenMP_CXX_FLAGS})

  find_package(AdaptiveCpp CONFIG REQUIRED)
  add_sycl_to_target(TARGET seissol-device-lib SOURCES ${DEVICE_SRC})

elseif("${DEVICE_BACKEND}" STREQUAL "oneapi")
  find_package(DpcppFlags REQUIRED)
  target_compile_options(seissol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-O3")
  target_compile_definitions(seissol-device-lib PRIVATE __DPCPP_COMPILER)
  target_link_libraries(seissol-device-lib PRIVATE dpcpp::device_flags ${GemmTools_LIBRARIES})
  target_link_libraries(seissol-common-properties INTERFACE dpcpp::interface)
endif()

if (USE_DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
target_compile_definitions(seissol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
endif()
