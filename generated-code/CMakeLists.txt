# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

set(GENERATED_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/")

add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${GENERATED_OUTPUT_DIR}")

set(GENERATED_FILES_FOR_SEISSOL ${GENERATED_OUTPUT_DIR}/subroutine.h
                                ${GENERATED_OUTPUT_DIR}/tensor.cpp
                                ${GENERATED_OUTPUT_DIR}/subroutine.cpp
                                ${GENERATED_OUTPUT_DIR}/gpulike_subroutine.cpp
                                ${GENERATED_OUTPUT_DIR}/tensor.h
                                ${GENERATED_OUTPUT_DIR}/init.cpp
                                ${GENERATED_OUTPUT_DIR}/init.h
                                ${GENERATED_OUTPUT_DIR}/kernel.h
                                ${GENERATED_OUTPUT_DIR}/kernel.cpp)

if (TESTING AND TESTING_GENERATED)
  set(GENERATED_FILES_FOR_SEISSOL ${GENERATED_FILES_FOR_SEISSOL} ${GENERATED_OUTPUT_DIR}/test-kernel.cpp)
endif()

file(MAKE_DIRECTORY ${GENERATED_OUTPUT_DIR})

if(PROXY_PYBINDING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

if(PREMULTIPLY_FLUX)
  set(PREMULTIPLY_FLUX_VALUE "--enable_premultiply_flux")
else()
  set(PREMULTIPLY_FLUX_VALUE "--disable_premultiply_flux")
endif()

find_package(GemmTools REQUIRED)

add_custom_command(
  COMMAND
  "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/generate.py"
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/matrices
     "--outputDir" ${GENERATED_OUTPUT_DIR}
     "--host_arch" ${HOST_ARCH_STR}
     "--device_codegen" ${DEVICE_CODEGEN}
     "--device_arch" ${DEVICE_ARCH_STR}
     "--device_backend" ${DEVICE_BACKEND}
     "--order" ${ORDER}
     "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" ${MEMORY_LAYOUT}
     "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
     "--gemm_tools" ${GEMM_TOOLS_LIST}
     "--drQuadRule" ${DR_QUAD_RULE}
     "--executable_libxsmm=${Libxsmm_executable_PROGRAM}"
     "--executable_pspamm=${PSpaMM_PROGRAM}"
     ${PREMULTIPLY_FLUX_VALUE} # boolean flag
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     DEPENDS
        build-time-make-directory
        ${CMAKE_CURRENT_SOURCE_DIR}/generate.py
        ${CMAKE_CURRENT_SOURCE_DIR}/kernels
        ${CMAKE_CURRENT_SOURCE_DIR}/matrices
        ${CMAKE_CURRENT_SOURCE_DIR}/config
    OUTPUT
        ${GENERATED_FILES_FOR_SEISSOL}
    COMMENT "Code generation for tensor operations..."
       )

add_custom_target(seissol-codegen ALL DEPENDS ${GENERATED_FILES_FOR_SEISSOL})

# Source code
add_library(seissol-kernel-lib
    ${GENERATED_OUTPUT_DIR}/kernel.cpp
    ${GENERATED_OUTPUT_DIR}/tensor.cpp
    ${GENERATED_OUTPUT_DIR}/subroutine.cpp
    ${GENERATED_OUTPUT_DIR}/init.cpp
)

if(GemmTools_INCLUDE_DIRS)
  target_include_directories(seissol-common-properties INTERFACE ${GemmTools_INCLUDE_DIRS})
endif()
if(GemmTools_LIBRARIES)
  target_link_libraries(seissol-common-properties INTERFACE ${GemmTools_LIBRARIES})
endif()
if(GemmTools_COMPILER_DEFINITIONS)
  target_compile_definitions(seissol-common-properties INTERFACE ${GemmTools_COMPILER_DEFINITIONS})
endif()

target_link_libraries(seissol-kernel-lib PUBLIC seissol-common-properties)

# Generated code does only work without red-zone.
if (HAS_REDZONE)
  set_source_files_properties(
    ${GENERATED_OUTPUT_DIR}/subroutine.cpp PROPERTIES COMPILE_FLAGS -mno-red-zone
  )
endif()
