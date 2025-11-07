# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff


add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode")

set(GENERATED_FILES_FOR_SEISSOL "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/subroutine.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/tensor.cpp"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/subroutine.cpp"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/gpulike_subroutine.cpp"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/tensor.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/init.cpp"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/init.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/kernel.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/kernel.cpp"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/test-kernel.cpp")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/)

if(PREMULTIPLY_FLUX)
  set(PREMULTIPLY_FLUX_VALUE "--enable_premultiply_flux")
else()
  set(PREMULTIPLY_FLUX_VALUE "--disable_premultiply_flux")
endif()

if (OVERRIDE_VECTORSIZE GREATER 0)
  set(OVERRIDE_VECTORSIZE_PARAM "--vectorsize=${OVERRIDE_VECTORSIZE}")
else()
  set(OVERRIDE_VECTORSIZE_PARAM "")
endif()

add_custom_command(
  COMMAND
  "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/generate.py"
     "--equations" ${EQUATIONS}
     "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/matrices
     "--outputDir" ${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode
     "--host_arch" ${HOST_ARCH}
     "--device_codegen" ${DEVICE_CODEGEN}
     "--device_arch" ${DEVICE_ARCH}
     "--device_backend" ${DEVICE_BACKEND}
     "--device_vendor" ${DEVICE_VENDOR}
     "--precision" ${PRECISION_PREFIX}
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
     ${OVERRIDE_VECTORSIZE_PARAM}
     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/codegen
     DEPENDS
        build-time-make-directory
        "${CMAKE_CURRENT_SOURCE_DIR}/generate.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/kernels"
        "${CMAKE_CURRENT_SOURCE_DIR}/matrices"
        "${CMAKE_CURRENT_SOURCE_DIR}/config"
    OUTPUT
        ${GENERATED_FILES_FOR_SEISSOL}
    COMMENT "Code generation for tensor operations..."
       )

add_custom_target(seissol-codegen ALL DEPENDS ${GENERATED_FILES_FOR_SEISSOL})
