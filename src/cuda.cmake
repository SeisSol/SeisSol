find_package(CUDA REQUIRED)

# Note: -std=c++14 because of cuda@10
set(CUDA_NVCC_FLAGS -std=c++14;
        -Xptxas -v;
        -arch=${DEVICE_ARCH};
        -DREAL_SIZE=${REAL_SIZE_IN_BYTES};
        --compiler-options ${EXTRA_CXX_FLAGS};
        -O3;)

set(DEVICE_SRC ${DEVICE_SRC}
               ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/cuda/PlasticityAux.cu)

set_source_files_properties(${DEVICE_SRC} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

cuda_add_library(SeisSol-device-lib STATIC ${DEVICE_SRC})
target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE}
                                                     ${CUDA_TOOLKIT_ROOT_DIR})
target_compile_options(SeisSol-device-lib PRIVATE ${EXTRA_CXX_FLAGS})
