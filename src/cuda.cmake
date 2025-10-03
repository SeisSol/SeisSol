# SPDX-FileCopyrightText: 2021 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

enable_language(CUDA)

set(DEVICE_SRC ${DEVICE_SRC}
    ${CMAKE_BINARY_DIR}/codegen/GeneratedCode/gpulike_subroutine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/cudahip/PlasticityAux.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/LinearCK/DeviceAux/cudahip/KernelsAux.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverCudaHip.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/PointSourceClusterCudaHip.cpp)

function(make_device_lib NAME FILES)

    add_library(${NAME} SHARED ${FILES})

    set_target_properties(${NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    set_source_files_properties(${FILES} PROPERTIES LANGUAGE CUDA)

    target_include_directories(${NAME} PUBLIC ${SEISSOL_DEVICE_INCLUDE})
    target_compile_features(${NAME} PRIVATE cxx_std_17)
    target_compile_features(${NAME} PRIVATE cuda_std_17)
    target_compile_options(${NAME} PRIVATE ${EXTRA_CXX_FLAGS})

    string(REPLACE "sm_" "" CUDA_DEVICE_ARCH "${DEVICE_ARCH}")
    set_target_properties(${NAME} PROPERTIES CUDA_ARCHITECTURES "${CUDA_DEVICE_ARCH}")

    target_compile_options(${NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
            -std=c++17;
            --expt-relaxed-constexpr;
            >)

    if (DEVICE_KERNEL_INFOPRINT)
        target_compile_options(${NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
            -Xptxas -v;
            >)
    endif()

    if (DEVICE_KERNEL_SAVETEMPS)
        target_compile_options(${NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
            --keep;
            >)
    endif()

    if (EXTRA_CXX_FLAGS)
        target_compile_options(${NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                --compiler-options ${EXTRA_CXX_FLAGS}
                >)
    endif()

endfunction()

make_device_lib(seissol-device-lib "${DEVICE_SRC}")
