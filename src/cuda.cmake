# SPDX-FileCopyrightText: 2021 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

enable_language(CUDA)

set(DEVICE_SRC ${DEVICE_SRC}
        ${CMAKE_BINARY_DIR}/generated-code/GeneratedCode/gpulike_subroutine.cpp
        Kernels/DeviceAux/cudahip/PlasticityAux.cpp
        Kernels/LinearCK/DeviceAux/cudahip/KernelsAux.cpp
        DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverCudaHip.cpp
        Kernels/PointSourceClusterCudaHip.cpp)

add_library(seissol-device-lib SHARED ${DEVICE_SRC})

set_target_properties(seissol-device-lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_source_files_properties(${DEVICE_SRC} PROPERTIES LANGUAGE CUDA)

target_link_libraries(seissol-device-lib PRIVATE seissol-common-properties)

target_include_directories(seissol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
target_compile_features(seissol-device-lib PRIVATE cxx_std_17)
target_compile_features(seissol-device-lib PRIVATE cuda_std_17)
target_compile_options(seissol-device-lib PRIVATE ${EXTRA_CXX_FLAGS})

if (USE_DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
target_compile_definitions(seissol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
endif()

string(REPLACE "sm_" "" CUDA_DEVICE_ARCH "${DEVICE_ARCH}")
set_target_properties(seissol-device-lib PROPERTIES CUDA_ARCHITECTURES "${CUDA_DEVICE_ARCH}")

target_compile_options(seissol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        -Xptxas -v;
        -std=c++17;
        --expt-relaxed-constexpr;
        >)
if (EXTRA_CXX_FLAGS)
    target_compile_options(seissol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
            --compiler-options ${EXTRA_CXX_FLAGS}
            >)
endif()

