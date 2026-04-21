# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

if(NOT Python3_EXECUTABLE)
  message(FATAL_ERROR "Python3_EXECUTABLE is required")
endif()

if(NOT SEISSOL_SOURCE_DIR)
  message(FATAL_ERROR "SEISSOL_SOURCE_DIR is required")
endif()

if(NOT SEISSOL_BINARY_DIR)
  message(FATAL_ERROR "SEISSOL_BINARY_DIR is required")
endif()

set(CodegenSourceDir "${SEISSOL_SOURCE_DIR}/codegen")
set(OutputDir "${SEISSOL_BINARY_DIR}/tests/GeneratedCodeMemoryLayoutFallback")
set(NonExistingAbsDense "/tmp/seissol-test-not-existing/dense.xml")

file(REMOVE_RECURSE "${OutputDir}")
file(MAKE_DIRECTORY "${OutputDir}")

execute_process(
  COMMAND
    "${Python3_EXECUTABLE}" "${CodegenSourceDir}/generate.py"
    "--equations" "elastic"
    "--matricesDir" "${CodegenSourceDir}/matrices"
    "--outputDir" "${OutputDir}"
    "--host_arch" "hsw"
    "--device_codegen" "none"
    "--device_arch" "none"
    "--device_backend" "none"
    "--device_vendor" "none"
    "--precision" "d"
    "--order" "6"
    "--numberOfMechanisms" "0"
    "--memLayout" "${NonExistingAbsDense}"
    "--multipleSimulations" "1"
    "--PlasticityMethod" "nb"
    "--gemm_tools" "none"
    "--drQuadRule" "stroud"
    "--disable_premultiply_flux"
  WORKING_DIRECTORY "${CodegenSourceDir}"
  RESULT_VARIABLE Result
  OUTPUT_VARIABLE Stdout
  ERROR_VARIABLE Stderr
)

if(NOT Result EQUAL 0)
  message(FATAL_ERROR
    "codegen/generate.py failed for absolute non-existing dense.xml path.\n"
    "stdout:\n${Stdout}\n"
    "stderr:\n${Stderr}")
endif()

string(FIND "${Stdout}" "Using the pre-defined memory layout config file dense.xml (cpu)" MessagePos)
if(MessagePos EQUAL -1)
  message(FATAL_ERROR
    "Expected fallback message was not found.\n"
    "stdout:\n${Stdout}\n"
    "stderr:\n${Stderr}")
endif()

file(GLOB GeneratedInitHeaders "${OutputDir}/equation-*/init.h")
if(GeneratedInitHeaders STREQUAL "")
  message(FATAL_ERROR
    "Expected generated init.h was not found under ${OutputDir}/equation-*/init.h")
endif()
