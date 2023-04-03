if ("${DEVICE_BACKEND}" STREQUAL "opensycl")

  set(DEVICE_SRC ${DEVICE_SRC}
          ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/sycl/PlasticityAux.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DeviceAux/sycl/KernelsAux.cpp)

  add_library(SeisSol-device-lib SHARED ${DEVICE_SRC})

  target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
  target_compile_definitions(SeisSol-device-lib PRIVATE DEVICE_HIPSYCL_LANG REAL_SIZE=${REAL_SIZE_IN_BYTES})
  target_link_libraries(SeisSol-device-lib PUBLIC ${Boost_LIBRARIES})

  find_package(OpenSYCLSettings REQUIRED)
  find_package(OpenSYCL REQUIRED)
  target_link_libraries(SeisSol-device-lib PRIVATE opensycl-settings::device_flags)
  add_sycl_to_target(TARGET SeisSol-device-lib SOURCES ${DEVICE_SRC})

elseif("${DEVICE_BACKEND}" STREQUAL "oneapi")
  set(DEVICE_SRC ${DEVICE_SRC}
                 ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/sycl/PlasticityAux.cpp)

  add_library(SeisSol-device-lib SHARED ${DEVICE_SRC})

  target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})

  target_compile_options(SeisSol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-O3")
  target_compile_definitions(SeisSol-device-lib PRIVATE DEVICE_ONEAPI_LANG REAL_SIZE=${REAL_SIZE_IN_BYTES})

  find_package(DpcppFlags REQUIRED)
  target_link_libraries(SeisSol-device-lib PRIVATE dpcpp::device_flags)
endif()
