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
    Core.cpp
    Cost.cpp
    RtcCpu.cpp
    RtcGpu.cpp
    RtcGpuDriver.cpp
    Validate.cpp
    Lower.cpp
    Branches.cpp
    Expand.cpp
    Interp.cpp
    Record.cpp
    Rewrite.cpp
    SderivFrontend.cpp
    SderivLibrary.cpp
    Sha256.cpp
    Types.cpp
)

target_link_libraries(seissol-expr PUBLIC seissol-common-properties)

# The interpreter is the reference the compiled backends are checked against bit for bit, and they
# are compiled without contraction of floating-point expressions; so is it. Otherwise a contraction
# a*b + c may become a fused multiply-add here and not there.
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
    target_compile_options(seissol-expr PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)
endif()
# icpx compiles with -fp-model=fast from -O2 on, i.e. -ffast-math with all its flags (no NaN, no
# infinities, no signed zero, reciprocals, contraction, approximate functions, reassociation): the
# reference would not be IEEE. The interpreter checks this at its first use (checkInterpreter()).
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(seissol-expr PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-fp-model=precise>)
endif()

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
#
# The driver calls the SYCL runtime from host code that is no device library, so it links the
# runtime itself: the executable has it only through the device libraries, and ld takes no symbol
# from a library that is no direct dependency ("DSO missing from command line").
if (WITH_GPU AND DEVICE_BACKEND STREQUAL "oneapi")
    get_filename_component(SEISSOL_CXX_COMPILER_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)
    find_library(SEISSOL_SYCL_LIBRARY sycl HINTS ${SEISSOL_CXX_COMPILER_DIR}/../lib)
    if (SEISSOL_SYCL_LIBRARY)
        target_link_libraries(seissol-expr PUBLIC ${SEISSOL_SYCL_LIBRARY})
        target_compile_definitions(seissol-expr PRIVATE SEISSOL_EXPR_HAVE_SYCL)
        message(STATUS "seissol-expr: OpenCL-through-SYCL backend enabled (${SEISSOL_SYCL_LIBRARY})")
    else()
        message(STATUS "seissol-expr: no SYCL runtime library; device programs fall back")
    endif()
endif()
