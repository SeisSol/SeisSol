# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff


add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode")

set(GENERATED_FILES_FOR_SEISSOL "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/subroutine.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/subroutine.cpp"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/gpulike_subroutine.cpp"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/tensor.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/init.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/kernel.h")

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

if (OVERRIDE_ALIGNMENT GREATER 0)
  set(OVERRIDE_ALIGNMENT_PARAM "--alignment=${OVERRIDE_ALIGNMENT}")
else()
  set(OVERRIDE_ALIGNMENT_PARAM "")
endif()

# common arguments between configure and build
set(GENERATE_ARGS "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/generate.py"
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
     "--numMechanisms" ${NUMBER_OF_MECHANISMS}
     "--memLayout" ${MEMORY_LAYOUT}
     "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
     "--PlasticityMethod" ${PLASTICITY_METHOD}
     "--gemm_tools" ${GEMM_TOOLS_LIST}
     "--drQuadRule" ${DR_QUAD_RULE}
     "--executable_libxsmm=${Libxsmm_executable_PROGRAM}"
     "--executable_pspamm=${PSpaMM_PROGRAM}"
     "--solver" ${SOLVER}
     ${PREMULTIPLY_FLUX_VALUE} # boolean flag
     ${OVERRIDE_VECTORSIZE_PARAM}
     ${OVERRIDE_ALIGNMENT_PARAM})

# find out which files will be generated
execute_process(COMMAND ${GENERATE_ARGS} "--mode" "collect"
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

file(READ "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/targets.json" SEISSOL_CODEGEN_TARGETS)

string(JSON TARGETCOUNT LENGTH ${SEISSOL_CODEGEN_TARGETS})
math(EXPR TARGETCOUNT "${TARGETCOUNT} - 1")

set(SEISSOL_CODEGEN_KERNELS)
set(SEISSOL_CODEGEN_TESTS)
set(SEISSOL_CODEGEN_HEADERS)

foreach(IDX RANGE ${TARGETCOUNT})
  string(JSON TARGETNAME MEMBER ${SEISSOL_CODEGEN_TARGETS} ${IDX})

  # into SEISSOL_CODEGEN_KERNELS, SEISSOL_CODEGEN_TESTS and SEISSOL_CODEGEN_HEADERS
  foreach(KIND kernels tests headers)
    string(TOUPPER ${KIND} KIND_LIST)
    string(JSON FILE_COUNT LENGTH ${SEISSOL_CODEGEN_TARGETS} "${TARGETNAME}" "${KIND}")
    # a target need not have files of every kind
    if (FILE_COUNT GREATER 0)
      math(EXPR LAST_FILE "${FILE_COUNT} - 1")
      foreach (IDX2 RANGE ${LAST_FILE})
        string(JSON TARGET_FILE GET ${SEISSOL_CODEGEN_TARGETS} "${TARGETNAME}" "${KIND}" ${IDX2})
        list(APPEND SEISSOL_CODEGEN_${KIND_LIST} "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/${TARGET_FILE}")
      endforeach()
    endif()
  endforeach()

endforeach()

# TODO: fix PARENT_SCOPE usage
set(SEISSOL_CODEGEN_TESTS ${SEISSOL_CODEGEN_TESTS} PARENT_SCOPE)

set(GENERATED_FILES_FOR_SEISSOL
      ${GENERATED_FILES_FOR_SEISSOL}
      ${SEISSOL_CODEGEN_KERNELS}
      ${SEISSOL_CODEGEN_TESTS}
      ${SEISSOL_CODEGEN_HEADERS})

add_custom_command(
  COMMAND
     ${GENERATE_ARGS} "--mode" "codegen"
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
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

add_library(seissol-kernel-lib
  ${SEISSOL_CODEGEN_KERNELS}
  ${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/subroutine.cpp
)
