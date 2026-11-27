# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Generated code does only work without red-zone.
if (HAS_REDZONE)
  set_source_files_properties(
      ${SEISSOL_CODEGEN_ROUTINES} PROPERTIES COMPILE_FLAGS -mno-red-zone
  )
endif()

target_compile_options(seissol-kernel-lib PRIVATE -fPIC)

if (SHARED)
  add_library(seissol-lib SHARED)
else()
  add_library(seissol-lib STATIC)
endif()

target_sources(seissol-lib PRIVATE
Simulation/SeisSol.cpp
)
target_link_libraries(seissol-lib PUBLIC seissol-config)

# include necessary kernel files (we can't include all of them right now, because of some undefined kernels + tensors)

# GPU code
if (WITH_GPU)
  # include cmake files will define seissol-device-lib target
  if ("${DEVICE_BACKEND}" STREQUAL "cuda" OR "${DEVICE_BACKEND}" STREQUAL "hip")
    set(DEVICE_SRC ${DEVICE_SRC}
      ${SEISSOL_CODEGEN_DEVICE}
      Core/Common/TypeId.cpp
      Core/Parallel/DeviceAux/cudahip/ClockAux.cpp
      Core/Parallel/DeviceAux/cudahip/GatherAux.cpp
      Core/Parallel/DeviceAux/cudahip/HostCopyAux.cpp
      Features/Plasticity/Solver/DeviceAux/cudahip/PlasticityAux.cpp
      WavePropagation/Solver/LinearCK/DeviceAux/cudahip/KernelsAux.cpp
      Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolverCudaHip.cpp
      Features/SourceTerm/Solver/PointSourceClusterCudaHip.cpp)
  elseif ("${DEVICE_BACKEND}" STREQUAL "hipsycl" OR "${DEVICE_BACKEND}" STREQUAL "acpp" OR "${DEVICE_BACKEND}" STREQUAL "oneapi")
    set(DEVICE_SRC ${DEVICE_SRC}
          ${SEISSOL_CODEGEN_DEVICE}
          Core/Common/TypeId.cpp
          Core/Parallel/DeviceAux/sycl/ClockAux.cpp
          Core/Parallel/DeviceAux/sycl/GatherAux.cpp
          Core/Parallel/DeviceAux/sycl/HostCopyAux.cpp
          Features/Plasticity/Solver/DeviceAux/sycl/PlasticityAux.cpp
          WavePropagation/Solver/LinearCK/DeviceAux/sycl/KernelsAux.cpp
          Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolverSycl.cpp
          Features/SourceTerm/Solver/PointSourceClusterSycl.cpp)
  endif()

  make_device_lib(seissol-device-lib "${DEVICE_SRC}")

  if (WARNINGS_AS_ERRORS)
    # mostly generated code, compiled by the GPU toolchains
    set_property(SOURCE ${DEVICE_SRC} APPEND PROPERTY COMPILE_OPTIONS
      $<$<COMPILE_LANGUAGE:CXX>:-Wno-error>)
  endif()

  target_link_libraries(seissol-device-lib PRIVATE seissol-common-properties)

  target_compile_options(seissol-device-lib PRIVATE -fPIC)
  target_include_directories(seissol-lib PRIVATE ${DEVICE_INCLUDE_DIRS})

  if (USE_DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
    target_compile_definitions(seissol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
  endif()
endif()
