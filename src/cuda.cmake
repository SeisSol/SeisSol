enable_language(CUDA)

set(DEVICE_SRC ${DEVICE_SRC}
        ${CMAKE_BINARY_DIR}/generated-code/generated-code/gpulike_subroutine.cpp
        ${CMAKE_SOURCE_DIR}/src/Kernels/DeviceAux/cuda/PlasticityAux.cu
        ${CMAKE_SOURCE_DIR}/src/Equations/elastic/Kernels/DeviceAux/cuda/KernelsAux.cu)

add_library(SeisSol-device-lib SHARED ${DEVICE_SRC})

set_target_properties(SeisSol-device-lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_source_files_properties(${DEVICE_SRC} PROPERTIES LANGUAGE CUDA)

target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
target_compile_features(SeisSol-device-lib PRIVATE cxx_std_17)
target_compile_features(SeisSol-device-lib PRIVATE cuda_std_17)
target_compile_options(SeisSol-device-lib PRIVATE ${EXTRA_CXX_FLAGS})
target_compile_definitions(SeisSol-device-lib PRIVATE ${HARDWARE_DEFINITIONS}
        CONVERGENCE_ORDER=${ORDER}
        NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
        NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
        ${DR_QUAD_RULE})

if (DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
target_compile_definitions(SeisSol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
endif()

string(REPLACE "sm_" "" CUDA_DEVICE_ARCH "${DEVICE_ARCH}")
set_target_properties(SeisSol-device-lib PROPERTIES CUDA_ARCHITECTURES "${CUDA_DEVICE_ARCH}")

target_compile_definitions(SeisSol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        -DREAL_SIZE=${REAL_SIZE_IN_BYTES}
        >)

target_compile_options(SeisSol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        -Xptxas -v;
        -std=c++17;
        --expt-relaxed-constexpr;
        >)
if (EXTRA_CXX_FLAGS)
    target_compile_options(SeisSol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
            --compiler-options ${EXTRA_CXX_FLAGS}
            >)
endif()
