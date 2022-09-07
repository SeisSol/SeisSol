enable_language(CUDA)
set(CMAKE_CUDA_STANDARD 14)

set(DEVICE_SRC ${DEVICE_SRC};
        ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp;
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/cuda/PlasticityAux.cu)

add_library(SeisSol-device-lib SHARED ${DEVICE_SRC})

set_target_properties(SeisSol-device-lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_source_files_properties(${DEVICE_SRC} PROPERTIES LANGUAGE CUDA)

target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
target_compile_features(SeisSol-device-lib PRIVATE cxx_std_14)
target_compile_options(SeisSol-device-lib PRIVATE ${EXTRA_CXX_FLAGS})

string(REPLACE "sm_" "" CUDA_DEVICE_ARCH "${DEVICE_ARCH}")
set_target_properties(SeisSol-device-lib PROPERTIES CUDA_ARCHITECTURES "${CUDA_DEVICE_ARCH}")

target_compile_definitions(SeisSol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        -DREAL_SIZE=${REAL_SIZE_IN_BYTES}
        >)

target_compile_options(SeisSol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        -Xptxas -v;
        --expt-relaxed-constexpr;
        >)
if (EXTRA_CXX_FLAGS)
    target_compile_options(SeisSol-device-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
            --compiler-options ${EXTRA_CXX_FLAGS}
            >)
endif()
