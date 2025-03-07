# SPDX-FileCopyrightText: 2021-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# ensure that we have set HIP_PATH
if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
        if (NOT DEFINED ENV{ROCM_PATH})
            # default location
            set(HIP_PATH "/opt/rocm" CACHE PATH "Path to which HIP has been installed")
        else()
            set(HIP_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which HIP has been installed")
        endif()
    else()
        set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
    endif()
endif()

# set the CMAKE_MODULE_PATH for the helper cmake files from HIP
set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" "${HIP_PATH}/lib/cmake/hip" ${CMAKE_MODULE_PATH})

find_package(HIP REQUIRED)

set(SEISSOL_HIPCC -DREAL_SIZE=${REAL_SIZE_IN_BYTES}; -std=c++17; -O3)

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
               ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DeviceAux/hip/KernelsAux.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolverCudaHip.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/PointSourceClusterCudaHip.cpp)


set_source_files_properties(${DEVICE_SRC} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)

hip_reset_flags()
hip_add_library(seissol-device-lib SHARED ${DEVICE_SRC}
        HIPCC_OPTIONS ${SEISSOL_HIPCC}
        NVCC_OPTIONS ${SEISSOL_NVCC})

target_link_libraries(seissol-device-lib PRIVATE seissol-common-properties)
target_include_directories(seissol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
set_property(TARGET seissol-device-lib PROPERTY HIP_ARCHITECTURES OFF)
target_compile_definitions(seissol-device-lib PRIVATE ${HARDWARE_DEFINITIONS}
        CONVERGENCE_ORDER=${ORDER}
        NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
        NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}
        ${DR_QUAD_RULE})
if (DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
target_compile_definitions(seissol-device-lib PRIVATE DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS)
endif()

if (IS_NVCC_PLATFORM)
    set_target_properties(seissol-device-lib PROPERTIES LINKER_LANGUAGE HIP)
    target_link_options(seissol-device-lib PRIVATE -arch=${DEVICE_ARCH})
else()
    target_link_libraries(seissol-device-lib PUBLIC ${HIP_PATH}/lib/libamdhip64.so)
endif()

