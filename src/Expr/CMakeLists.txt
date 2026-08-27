# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-expr OBJECT
    Ir.cpp
    Program.cpp
    Binding.cpp
    Backend.cpp
    Codegen.cpp
    Cost.cpp
    RtcCpu.cpp
    RtcGpu.cpp
    RtcGpuDriver.cpp
    Validate.cpp
    Lower.cpp
    Interp.cpp
    SderivFrontend.cpp
)

target_link_libraries(seissol-expr PUBLIC seissol-common-properties)

# The compiled CPU backend loads its artifact with dlopen. It degrades to the
# interpreter at run time when the toolchain is unusable, so this is the only
# build-time requirement it has.
target_link_libraries(seissol-expr PUBLIC ${CMAKE_DL_LIBS})
target_link_libraries(seissol-lib PUBLIC seissol-expr)

# Runtime compilation for the device backends.
#
# Optional in both directions: without a runtime compiler the driver layer
# compiles to a stub that returns nullptr, and makeKernel falls back with a
# warning naming the reason. The code generator itself needs no toolkit and is
# tested without one -- the emitted source is compiled for the host and checked
# against the interpreter.
if (WITH_GPU AND DEVICE_BACKEND STREQUAL "cuda")
    find_package(CUDAToolkit)
    if (TARGET CUDA::nvrtc AND TARGET CUDA::cuda_driver)
        target_link_libraries(seissol-expr PUBLIC CUDA::nvrtc CUDA::cuda_driver)
        target_compile_definitions(seissol-expr PRIVATE SEISSOL_EXPR_HAVE_NVRTC)
        message(STATUS "seissol-expr: NVRTC backend enabled")
    else()
        message(STATUS "seissol-expr: no NVRTC; device programs fall back")
    endif()
endif()

if (WITH_GPU AND DEVICE_BACKEND STREQUAL "hip")
    find_package(hip QUIET)
    find_library(SEISSOL_HIPRTC_LIBRARY hiprtc HINTS ${HIP_ROOT_DIR}/lib $ENV{ROCM_PATH}/lib)
    if (SEISSOL_HIPRTC_LIBRARY)
        target_link_libraries(seissol-expr PUBLIC ${SEISSOL_HIPRTC_LIBRARY})
        target_compile_definitions(seissol-expr PRIVATE SEISSOL_EXPR_HAVE_HIPRTC)
        message(STATUS "seissol-expr: HIPRTC backend enabled")
    else()
        message(STATUS "seissol-expr: no HIPRTC; device programs fall back")
    endif()
endif()

# The Intel path: OpenCL C source enqueued through SYCL's kernel_compiler
# extension. Whether a given DEVICE accepts OpenCL C source is queried at run
# time, so this only says the compiler and headers are present.
if (WITH_GPU AND DEVICE_BACKEND STREQUAL "oneapi")
    target_compile_definitions(seissol-expr PRIVATE SEISSOL_EXPR_HAVE_SYCL)
    message(STATUS "seissol-expr: OpenCL-through-SYCL backend enabled")
endif()
