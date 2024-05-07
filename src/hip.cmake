
# ensure that we have an set HIP_PATH
if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
        set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
    else()
        set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
    endif()
endif()

# set the CMAKE_MODULE_PATH for the helper cmake files from HIP
set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" "${HIP_PATH}/lib/cmake/hip" ${CMAKE_MODULE_PATH})

find_package(HIP REQUIRED)

# Note: -std=c++14 because of cuda@10
set(SEISSOL_HIPCC -DREAL_SIZE=${REAL_SIZE_IN_BYTES}; -std=c++14; -O3)

set(IS_NVCC_PLATFORM OFF)
if (DEFINED ENV{HIP_PLATFORM})
    if ($ENV{HIP_PLATFORM} STREQUAL "nvidia")
        set(IS_NVCC_PLATFORM ON)
    endif()
endif()

if (IS_NVCC_PLATFORM)
   set(SEISSOL_NVCC -arch=${DEVICE_ARCH};
                    -dc;
                    --expt-relaxed-constexpr;
                    --compiler-options -fPIC;
                    -DCUDA_UNDERHOOD)
else()
    set(SEISSOL_HIPCC ${SEISSOL_HIPCC} --offload-arch=${DEVICE_ARCH})
endif()


set(CMAKE_HIP_CREATE_SHARED_LIBRARY
"${HIP_HIPCC_CMAKE_LINKER_HELPER} \
${HCC_PATH} \
<CMAKE_SHARED_LIBRARY_CXX_FLAGS> \
<LANGUAGE_COMPILE_FLAGS> \
<LINK_FLAGS> \
<CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS> -o <TARGET> \
<OBJECTS> \
<LINK_LIBRARIES>")

set(DEVICE_SRC ${DEVICE_SRC}
               ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/hip/PlasticityAux.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DeviceAux/hip/KernelsAux.cpp)


set_source_files_properties(${DEVICE_SRC} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)

hip_reset_flags()
hip_add_library(SeisSol-device-lib SHARED ${DEVICE_SRC}
        HIPCC_OPTIONS ${SEISSOL_HIPCC}
        NVCC_OPTIONS ${SEISSOL_NVCC})

target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
set_property(TARGET SeisSol-device-lib PROPERTY HIP_ARCHITECTURES OFF)


if (IS_NVCC_PLATFORM)
    set_target_properties(SeisSol-device-lib PROPERTIES LINKER_LANGUAGE HIP)
    target_link_options(SeisSol-device-lib PRIVATE -arch=${DEVICE_ARCH})
else()
    target_link_libraries(SeisSol-device-lib PUBLIC ${HIP_PATH}/lib/libamdhip64.so)
endif()
