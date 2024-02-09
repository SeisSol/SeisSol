set(DEVICE_SRC ${DEVICE_SRC}
        ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/omptarget/PlasticityAux.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DeviceAux/omptarget/KernelsAux.cpp)

        add_library(SeisSol-device-lib SHARED ${DEVICE_SRC})
        target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
        target_compile_definitions(SeisSol-device-lib PRIVATE DEVICE_HIPSYCL_LANG REAL_SIZE=${REAL_SIZE_IN_BYTES})
        find_package(OpenMP REQUIRED)
        target_compile_options(SeisSol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-fPIC")

        find_package(OmpTargetFlags REQUIRED)
        target_link_libraries(SeisSol-device-lib PRIVATE omptarget::device_flags)
