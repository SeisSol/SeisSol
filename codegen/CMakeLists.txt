# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff


add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode")

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

# the unit tests of the kernels, only where they are compiled
if (TESTING_GENERATED)
  set(UNIT_TESTS_VALUE "doctest")
else()
  set(UNIT_TESTS_VALUE "none")
endif()

# The routines of the device are spread over several files, which compile at the same time.
if (NOT WITH_GPU)
  set(DEVICE_SHARDS 1)
elseif (DEVICE_CODEGEN_SHARDS GREATER 0)
  set(DEVICE_SHARDS ${DEVICE_CODEGEN_SHARDS})
else()
  math(EXPR DEVICE_SHARDS "4 * (${SEISSOL_CONFIG_LAST} + 1)")
endif()

# common arguments between configure and build
# the configurations to generate (cmake/configs.cmake), each under the C++ type it is generated for
set(SEISSOL_CONFIGS_JSON "[")
foreach(IDX RANGE ${SEISSOL_CONFIG_LAST})
  foreach(FIELD TYPES MATERIALS SOLVERS MECHANISMS ORDERS PRECISIONS DRQUADRULES SIMULATIONS)
    list(GET SEISSOL_CONFIG_${FIELD} ${IDX} CONFIG_${FIELD})
  endforeach()
  string(SUBSTRING ${CONFIG_PRECISIONS} 0 1 CONFIG_PRECISION_PREFIX)
  if (IDX GREATER 0)
    string(APPEND SEISSOL_CONFIGS_JSON ",")
  endif()
  string(APPEND SEISSOL_CONFIGS_JSON
    "\n  {\"key\": \"seissol::${CONFIG_TYPES}\", \"equations\": \"${CONFIG_MATERIALS}\", "
    "\"solver\": \"${CONFIG_SOLVERS}\", \"numMechanisms\": ${CONFIG_MECHANISMS}, "
    "\"order\": ${CONFIG_ORDERS}, \"precision\": \"${CONFIG_PRECISION_PREFIX}\", "
    "\"drQuadRule\": \"${CONFIG_DRQUADRULES}\", "
    "\"multipleSimulations\": ${CONFIG_SIMULATIONS}}")
endforeach()
string(APPEND SEISSOL_CONFIGS_JSON "\n]\n")
# written only when its content changes: the code generation depends on it, and a rerun of CMake
# should not regenerate and rebuild everything
file(CONFIGURE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/configs.json"
     CONTENT "${SEISSOL_CONFIGS_JSON}" @ONLY)

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
     "--configs" "${CMAKE_CURRENT_BINARY_DIR}/configs.json"
     ${PREMULTIPLY_FLUX_VALUE} # boolean flag
     ${OVERRIDE_VECTORSIZE_PARAM}
     ${OVERRIDE_ALIGNMENT_PARAM}
     "--unit_tests" ${UNIT_TESTS_VALUE}
     "--device_shards" ${DEVICE_SHARDS})

# find out which files will be generated
execute_process(COMMAND ${GENERATE_ARGS} "--mode" "collect"
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

file(READ "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/targets.json" SEISSOL_CODEGEN_TARGETS)

string(JSON TARGETCOUNT LENGTH ${SEISSOL_CODEGEN_TARGETS})
math(EXPR TARGETCOUNT "${TARGETCOUNT} - 1")

set(SEISSOL_CODEGEN_KERNELS)
set(SEISSOL_CODEGEN_DEVICE)
set(SEISSOL_CODEGEN_ROUTINES)
set(SEISSOL_CODEGEN_TESTS)
set(SEISSOL_CODEGEN_HEADERS)

foreach(IDX RANGE ${TARGETCOUNT})
  string(JSON TARGETNAME MEMBER ${SEISSOL_CODEGEN_TARGETS} ${IDX})

  # into SEISSOL_CODEGEN_KERNELS, SEISSOL_CODEGEN_DEVICE, SEISSOL_CODEGEN_TESTS and
  # SEISSOL_CODEGEN_HEADERS
  foreach(KIND kernels device tests headers)
    string(TOUPPER ${KIND} KIND_LIST)
    # a target need not have files of every kind, nor list every kind
    string(JSON FILE_COUNT ERROR_VARIABLE NO_KIND LENGTH ${SEISSOL_CODEGEN_TARGETS} "${TARGETNAME}" "${KIND}")
    if (NOT NO_KIND AND FILE_COUNT GREATER 0)
      math(EXPR LAST_FILE "${FILE_COUNT} - 1")
      foreach (IDX2 RANGE ${LAST_FILE})
        string(JSON TARGET_FILE GET ${SEISSOL_CODEGEN_TARGETS} "${TARGETNAME}" "${KIND}" ${IDX2})
        list(APPEND SEISSOL_CODEGEN_${KIND_LIST} "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/${TARGET_FILE}")
        if (TARGETNAME STREQUAL "routines" AND KIND STREQUAL "kernels")
          list(APPEND SEISSOL_CODEGEN_ROUTINES "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/${TARGET_FILE}")
        endif()
      endforeach()
    endif()
  endforeach()

endforeach()

# TODO: fix PARENT_SCOPE usage
set(SEISSOL_CODEGEN_TESTS ${SEISSOL_CODEGEN_TESTS} PARENT_SCOPE)
set(SEISSOL_CODEGEN_ROUTINES ${SEISSOL_CODEGEN_ROUTINES} PARENT_SCOPE)
set(SEISSOL_CODEGEN_DEVICE ${SEISSOL_CODEGEN_DEVICE} PARENT_SCOPE)

# Every configuration is generated by a command of its own, which the build runs at the same time
# as the others; the shared step reads what they hand to it (see STEP_ALL in generate.py).
file(READ "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/steps.json" SEISSOL_CODEGEN_STEPS)
string(JSON STEPCOUNT LENGTH ${SEISSOL_CODEGEN_STEPS})
math(EXPR STEPCOUNT "${STEPCOUNT} - 1")

set(GENERATED_FILES_FOR_SEISSOL)
foreach(IDX RANGE ${STEPCOUNT})
  string(JSON STEP_TARGET GET ${SEISSOL_CODEGEN_STEPS} ${IDX} "target")
  foreach(FIELD outputs inputs)
    set(STEP_${FIELD})
    string(JSON FILE_COUNT LENGTH ${SEISSOL_CODEGEN_STEPS} ${IDX} "${FIELD}")
    if (FILE_COUNT GREATER 0)
      math(EXPR LAST_FILE "${FILE_COUNT} - 1")
      foreach (IDX2 RANGE ${LAST_FILE})
        string(JSON STEP_FILE GET ${SEISSOL_CODEGEN_STEPS} ${IDX} "${FIELD}" ${IDX2})
        list(APPEND STEP_${FIELD} "${CMAKE_CURRENT_BINARY_DIR}/GeneratedCode/${STEP_FILE}")
      endforeach()
    endif()
  endforeach()

  if (STEP_TARGET STREQUAL "__shared__")
    set(STEP_COMMENT "Code generation for what the configurations share")
  else()
    set(STEP_COMMENT "Code generation for ${STEP_TARGET}")
  endif()

  add_custom_command(
    COMMAND
       ${GENERATE_ARGS} "--mode" "codegen" "--codegen_target" "${STEP_TARGET}"
       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
       DEPENDS
          build-time-make-directory
          "${CMAKE_CURRENT_SOURCE_DIR}/generate.py"
          "${CMAKE_CURRENT_SOURCE_DIR}/kernels"
          "${CMAKE_CURRENT_SOURCE_DIR}/matrices"
          "${CMAKE_CURRENT_SOURCE_DIR}/config"
          "${CMAKE_CURRENT_BINARY_DIR}/configs.json"
          ${STEP_inputs}
      OUTPUT
          ${STEP_outputs}
      COMMENT "${STEP_COMMENT}"
         )
  list(APPEND GENERATED_FILES_FOR_SEISSOL ${STEP_outputs})
endforeach()

add_custom_target(seissol-codegen ALL DEPENDS ${GENERATED_FILES_FOR_SEISSOL})

add_library(seissol-kernel-lib
  ${SEISSOL_CODEGEN_KERNELS}
)
