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
SeisSol.cpp
)
target_link_libraries(seissol-lib PUBLIC seissol-config)

# include necessary kernel files (we can't include all of them right now, because of some undefined kernels + tensors)

# the kernels of the solvers the configurations advance their cells with (cmake/configs.cmake);
# stp takes the local and neighbor kernels of linearck
if ("linearck" IN_LIST SEISSOL_SOLVERS_USED OR "stp" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCK/Local.cpp
    Kernels/LinearCK/Neighbor.cpp
    )
endif()
if ("linearck" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCK/Time.cpp
    )
  target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_LINEARCK)
endif()
if ("linearckanelastic" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-lib PRIVATE
    Kernels/LinearCKAnelastic/Neighbor.cpp
    Kernels/LinearCKAnelastic/Local.cpp
    Kernels/LinearCKAnelastic/Time.cpp
    )
  target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_LINEARCKANELASTIC)
endif()
if ("stp" IN_LIST SEISSOL_SOLVERS_USED)
  target_sources(seissol-lib PRIVATE
    Kernels/STP/Time.cpp
    )
  target_compile_definitions(seissol-common-properties INTERFACE SEISSOL_KERNELS_STP)
endif()

# the material headers of the equation of the first configuration
if ("${EQUATIONS}" STREQUAL "elastic" OR "${EQUATIONS}" STREQUAL "acoustic" OR "${EQUATIONS}" STREQUAL "anisotropic")
  target_include_directories(seissol-common-properties INTERFACE Equations/elastic)
elseif ("${EQUATIONS}" STREQUAL "viscoelastic" OR "${EQUATIONS}" STREQUAL "viscoacoustic")
  target_include_directories(seissol-common-properties INTERFACE Equations/viscoelastic)
elseif ("${EQUATIONS}" STREQUAL "poroelastic")
  target_include_directories(seissol-common-properties INTERFACE Equations/poroelastic)
endif()


# GPU code
if (WITH_GPU)
  # include cmake files will define seissol-device-lib target
  if ("${DEVICE_BACKEND}" STREQUAL "cuda" OR "${DEVICE_BACKEND}" STREQUAL "hip")
    set(DEVICE_SRC ${DEVICE_SRC}
      ${SEISSOL_CODEGEN_DEVICE}
      Kernels/DeviceAux/cudahip/PlasticityAux.cpp
      Kernels/LinearCK/DeviceAux/cudahip/KernelsAux.cpp
      DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverCudaHip.cpp
      Kernels/PointSourceClusterCudaHip.cpp)
  elseif ("${DEVICE_BACKEND}" STREQUAL "hipsycl" OR "${DEVICE_BACKEND}" STREQUAL "acpp" OR "${DEVICE_BACKEND}" STREQUAL "oneapi")
    set(DEVICE_SRC ${DEVICE_SRC}
          ${SEISSOL_CODEGEN_DEVICE}
          Kernels/DeviceAux/sycl/PlasticityAux.cpp
          Kernels/LinearCK/DeviceAux/sycl/KernelsAux.cpp
          DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverSycl.cpp
          Kernels/PointSourceClusterSycl.cpp)
  endif()

  make_device_lib(seissol-device-lib "${DEVICE_SRC}")

  target_link_libraries(seissol-device-lib PRIVATE seissol-common-properties)

  target_compile_options(seissol-device-lib PRIVATE -fPIC)
  target_include_directories(seissol-lib PRIVATE ${DEVICE_INCLUDE_DIRS})

  if (USE_DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
    target_compile_definitions(seissol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
  endif()
endif()
