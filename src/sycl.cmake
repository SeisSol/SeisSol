if ("${DEVICE_BACKEND}" STREQUAL "hipsycl")

  set(DEVICE_SRC ${DEVICE_SRC}
          ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/sycl/PlasticityAux.cpp)

  add_library(SeisSol-device-lib SHARED ${DEVICE_SRC})

  find_package(Boost REQUIRED COMPONENTS context fiber)
  if (DEVICE_ARCH MATCHES "sm_*")
    find_package(CUDA REQUIRED)
    set(HIPSYCL_TARGETS "cuda:${DEVICE_ARCH}")
    target_include_directories(SeisSol-device-lib PRIVATE ${CUDA_TOOLKIT_ROOT_DIR})
    target_link_libraries(SeisSol-device-lib PRIVATE -lcuda ${CUDA_LIBRARIES})
  elseif(DEVICE_ARCH MATCHES "gfx*")
    set(HIP_COMPILER hcc)
    find_package(HIP REQUIRED)
    set(HIPSYCL_TARGETS "hip:${DEVICE_ARCH}")
  else()
    set(HIPSYCL_TARGETS "${DEVICE_BACKEND}:${DEVICE_ARCH}")
  endif()

  target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})
  target_compile_definitions(SeisSol-device-lib PRIVATE DEVICE_HIPSYCL_LANG REAL_SIZE=${REAL_SIZE_IN_BYTES})
  target_link_libraries(SeisSol-device-lib PUBLIC ${Boost_LIBRARIES})

  find_package(OpenMP REQUIRED)
  target_compile_options(SeisSol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-fPIC" ${OpenMP_CXX_FLAGS})

  find_package(hipSYCL CONFIG REQUIRED)
  add_sycl_to_target(TARGET SeisSol-device-lib SOURCES ${DEVICE_SRC})

elseif("${DEVICE_BACKEND}" STREQUAL "oneapi")
  set(DEVICE_SRC ${DEVICE_SRC}
                 ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/sycl/PlasticityAux.cpp)

  add_library(SeisSol-device-lib STATIC ${DEVICE_SRC})

  target_include_directories(SeisSol-device-lib PUBLIC ${SEISSOL_DEVICE_INCLUDE})

  target_compile_options(SeisSol-device-lib PRIVATE ${EXTRA_CXX_FLAGS} "-O3")
  target_compile_definitions(SeisSol-device-lib PRIVATE DEVICE_ONEAPI_LANG REAL_SIZE=${REAL_SIZE_IN_BYTES})

  if ("$ENV{PREFERRED_DEVICE_TYPE}" STREQUAL "FPGA")
    message(NOTICE "FPGA is used as target device, compilation will take several hours to complete!")
    target_compile_options(SeisSol-device-lib PRIVATE -fsycl -fintelfpga -fsycl-unnamed-lambda)
    set_target_properties(SeisSol-device-lib PROPERTIES LINK_FLAGS -fsycl -fintelfpga -Xshardware)
  elseif("$ENV{PREFERRED_DEVICE_TYPE}" STREQUAL "GPU")
    if (${DEVICE_ARCH} MATCHES "sm_*")
      set(DEVICE_FLAGS -std=c++17
              -fsycl
              -Xsycl-target-backend
              --cuda-gpu-arch=${DEVICE_ARCH}
              -fsycl-targets=nvptx64-nvidia-cuda)
      target_compile_options(SeisSol-device-lib PRIVATE ${DEVICE_FLAGS})
      target_link_libraries(SeisSol-device-lib PRIVATE ${DEVICE_FLAGS})
    else()
      target_compile_options(SeisSol-device-lib PRIVATE "-fsycl" "-fsycl-targets=spir64_gen-unknown-unknown-sycldevice" "-fsycl-unnamed-lambda")
      set_target_properties(SeisSol-device-lib PROPERTIES LINK_FLAGS "-fsycl -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xs \"-device ${DEVICE_ARCH}\"")

      target_link_libraries(SeisSol-device-lib PUBLIC sycl "-fsycl" "-fsycl -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xs \"-device ${DEVICE_ARCH}\"")
      target_link_libraries(SeisSol-lib PUBLIC SeisSol-device-lib "-fsycl -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xs \"-device ${DEVICE_ARCH}\"")
      endif()
  elseif ("$ENV{PREFERRED_DEVICE_TYPE}" STREQUAL "CPU")
    target_compile_options(SeisSol-device-lib PRIVATE "-fsycl" "-fsycl-targets=spir64_x86_64-unknown-unknown-sycldevice" "-fsycl-unnamed-lambda")
    set_target_properties(SeisSol-device-lib PROPERTIES LINK_FLAGS "-fsycl -fsycl-targets=spir64_x86_64-unknown-unknown-sycldevice -Xs \"-march=${DEVICE_ARCH}\"")

    target_link_libraries(SeisSol-device-lib PUBLIC sycl "-fsycl" "-fsycl -fsycl-targets=spir64_x86_64-unknown-unknown-sycldevice -Xs \"-march ${DEVICE_ARCH}\"")
    target_link_libraries(SeisSol-lib PUBLIC SeisSol-device-lib "-fsycl -sycl-targets=spir64_x86_64-unknown-unknown-sycldevice -Xs \"-march ${DEVICE_ARCH}\"")
  else()
    message(FATAL_ERROR "please set PREFERRED_DEVICE_TYPE type to GPU, FPGA, or CPU in order to activate AOT "
            "compilation. If AOT is not performed, unnamed lambdas will cause errors at runtime")
  endif()
endif()
