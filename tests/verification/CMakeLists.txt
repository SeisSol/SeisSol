# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# ---------------------------------------------------------------------------
# Verification tests on generated mini meshes.
#
# The case list below is the declared matrix. Every build registers all of it
# and skips what it cannot run, so that a configuration missing from a build
# shows up as a skip instead of disappearing. The coverage test at the end
# turns those skips into a statement about what has been covered overall.
# ---------------------------------------------------------------------------

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(SEISSOL_MINITEST_TIMEOUT 300 CACHE STRING
    "Timeout in seconds for a single verification case")
set(SEISSOL_MINITEST_MPIEXEC_FLAGS "" CACHE STRING
    "Extra launcher flags for the verification cases, e.g. --oversubscribe on a \
machine with fewer cores than the rank counts the cases ask for")
set(SEISSOL_VERIFY_SCRIPT "" CACHE FILEPATH
    "Path to verify.py; probed under the source tree when empty")

if (NOT SEISSOL_VERIFY_SCRIPT)
  foreach (candidate scripts/verify.py .ci/verify.py tools/verify.py verify.py)
    if (EXISTS "${PROJECT_SOURCE_DIR}/${candidate}")
      set(SEISSOL_VERIFY_SCRIPT "${PROJECT_SOURCE_DIR}/${candidate}"
          CACHE FILEPATH "" FORCE)
      break()
    endif()
  endforeach()
endif()
if (SEISSOL_VERIFY_SCRIPT)
  message(STATUS "Verification: comparison script at ${SEISSOL_VERIFY_SCRIPT}")
else()
  message(STATUS "Verification: no verify.py found, comparison cases will skip")
endif()

# ---------------------------------------------------------------------------
# What this build can do.
#
# Only these axes are fixed at compile time. Plasticity, the friction law, the
# boundary conditions, the flux and the time stepping are all runtime
# parameters, so one build covers most of the matrix and the skip mechanism
# has to carry less weight than the number of axes suggests.
# ---------------------------------------------------------------------------

foreach (flag MPI OPENMP)
  if (${flag})
    set(SEISSOL_CAP_${flag} "true")
  else()
    set(SEISSOL_CAP_${flag} "false")
  endif()
endforeach()
if (DEVICE_BACKEND)
  set(SEISSOL_CAP_DEVICE "${DEVICE_BACKEND}")
else()
  set(SEISSOL_CAP_DEVICE "none")
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Capabilities.json.in
    ${CMAKE_CURRENT_BINARY_DIR}/capabilities.json
    @ONLY
)

# ---------------------------------------------------------------------------
# Meshes. Generated rather than checked in, as a setup fixture so that the
# cases stay independent and can run in parallel.
# ---------------------------------------------------------------------------

set(SEISSOL_MINIMESH_DIR ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(SEISSOL_MINIMESHES
    box6 box48 skew48 fault48 faultskew48 mat48 periodicx72 periodic162 open96)

add_test(NAME verification.minimesh.self-test
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/minimesh.py self-test)
set_tests_properties(verification.minimesh.self-test PROPERTIES
                     LABELS "verification;smoke")

foreach (mesh ${SEISSOL_MINIMESHES})
  add_test(NAME verification.mesh.${mesh}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/minimesh.py generate ${mesh}
                   -o ${SEISSOL_MINIMESH_DIR}/${mesh}.puml.h5
                   --manifest ${SEISSOL_MINIMESH_DIR}/${mesh}.json)
  set_tests_properties(verification.mesh.${mesh} PROPERTIES
                       FIXTURES_SETUP minimesh.${mesh}
                       LABELS "verification;mesh")
endforeach()

# ---------------------------------------------------------------------------
# seissol_add_minitest(NAME <name> MESH <mesh> PARAMETERS <file>
#                      [RANKS n] [THREADS n] [STEPS n] [TIMESTEP dt]
#                      [REQUIRE key=value ...] [SET Section.Key=Value ...]
#                      [FILES name ...] [ACTIVITY column:threshold ...]
#                      [LABELS label ...])
#
# REQUIRE is matched against capabilities.json by the runner, not here: a case
# that this build cannot run has to remain a visible skip.
# ---------------------------------------------------------------------------

define_property(GLOBAL PROPERTY SEISSOL_DECLARED_CASES
                BRIEF_DOCS "declared verification cases"
                FULL_DOCS "declared verification cases")

function(seissol_add_minitest)
  cmake_parse_arguments(MT
      ""
      "NAME;MESH;PARAMETERS;RANKS;THREADS;STEPS;TIMESTEP;THRESHOLDS"
      "REQUIRE;SET;FILES;ACTIVITY;LABELS"
      ${ARGN})
  if (NOT MT_NAME OR NOT MT_MESH OR NOT MT_PARAMETERS)
    message(FATAL_ERROR "seissol_add_minitest needs NAME, MESH and PARAMETERS")
  endif()
  if (NOT MT_MESH IN_LIST SEISSOL_MINIMESHES)
    message(FATAL_ERROR "unknown mini mesh '${MT_MESH}' in case '${MT_NAME}'")
  endif()
  if (NOT MT_RANKS)
    set(MT_RANKS 1)
  endif()
  if (NOT MT_THREADS)
    set(MT_THREADS 1)
  endif()

  # a case that asks for several ranks or threads requires a build that has
  # them; deriving it here keeps the declaration from having to repeat itself,
  # and makes such a case a skip rather than a launcher error
  if (MT_RANKS GREATER 1)
    list(APPEND MT_REQUIRE "mpi=true")
  endif()
  if (MT_THREADS GREATER 1)
    list(APPEND MT_REQUIRE "openmp=true")
  endif()
  list(REMOVE_DUPLICATES MT_REQUIRE)
  list(REMOVE_DUPLICATES MT_FILES)

  set(arguments
      run
      --name ${MT_NAME}
      --binary $<TARGET_FILE:seissol-bin>
      --capabilities ${CMAKE_CURRENT_BINARY_DIR}/capabilities.json
      --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/cases
      --work-dir ${CMAKE_CURRENT_BINARY_DIR}/work/${MT_NAME}
      --result ${CMAKE_CURRENT_BINARY_DIR}/results/${MT_NAME}.json
      --mesh ${SEISSOL_MINIMESH_DIR}/${MT_MESH}.puml.h5
      --parameters ${MT_PARAMETERS}
      --ranks ${MT_RANKS}
      --threads ${MT_THREADS})
  if (MT_STEPS)
    list(APPEND arguments --steps ${MT_STEPS})
  endif()
  if (MT_TIMESTEP)
    list(APPEND arguments --timestep ${MT_TIMESTEP})
  endif()
  foreach (entry ${MT_SET})
    list(APPEND arguments --set ${entry})
  endforeach()
  foreach (entry ${MT_FILES})
    list(APPEND arguments --file ${entry})
  endforeach()
  foreach (entry ${MT_ACTIVITY})
    list(APPEND arguments --activity ${entry})
  endforeach()
  foreach (entry ${MT_REQUIRE})
    list(APPEND arguments --require ${entry})
  endforeach()

  if (MT_RANKS GREATER 1 AND MPIEXEC_EXECUTABLE)
    # the equals form is required: these values start with a dash and would
    # otherwise be taken for the next option
    list(APPEND arguments
         "--mpiexec=${MPIEXEC_EXECUTABLE}"
         "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}")
    foreach (flag ${MPIEXEC_PREFLAGS} ${SEISSOL_MINITEST_MPIEXEC_FLAGS})
      list(APPEND arguments "--mpiexec-flag=${flag}")
    endforeach()
  endif()

  if (MT_THRESHOLDS)
    list(APPEND arguments --thresholds ${CMAKE_CURRENT_SOURCE_DIR}/cases/${MT_THRESHOLDS})
  endif()

  add_test(NAME verification.${MT_NAME}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py ${arguments})

  math(EXPR slots "${MT_RANKS} * ${MT_THREADS}")
  set_tests_properties(verification.${MT_NAME} PROPERTIES
      # the runner reports a configuration this build lacks as skipped, so that
      # it is neither a failure nor silently absent
      SKIP_RETURN_CODE 77
      FIXTURES_SETUP case.${MT_NAME}
      FIXTURES_REQUIRED "minimesh.${MT_MESH};verification"
      PROCESSORS ${slots}
      TIMEOUT ${SEISSOL_MINITEST_TIMEOUT}
      LABELS "verification;${MT_LABELS}")

  string(REPLACE ";" "\", \"" requires_json "${MT_REQUIRE}")
  if (requires_json)
    set(requires_json "\"${requires_json}\"")
  endif()
  set_property(GLOBAL APPEND_STRING PROPERTY SEISSOL_DECLARED_CASES
      "    {\"name\": \"${MT_NAME}\", \"mesh\": \"${MT_MESH}\", \"requires\": [${requires_json}]},\n")
endfunction()

# ---------------------------------------------------------------------------
# Stage 0: does it run at all, and did it compute anything.
#
# A fixed time step width and a step count keep the runs comparable between
# builds; without both, the step size follows the mesh and the material and
# nothing can be held against anything else.
# ---------------------------------------------------------------------------

set(SEISSOL_WAVE_MATERIALS "equations=elastic|anisotropic|viscoelastic2|poroelastic|acoustic")

foreach (mesh box48 skew48)
  seissol_add_minitest(
      NAME smoke-${mesh}
      MESH ${mesh}
      PARAMETERS wave.par
      FILES material.yaml receiver.dat
      STEPS 20
      TIMESTEP 1e-4
      REQUIRE ${SEISSOL_WAVE_MATERIALS}
      ACTIVITY "v1:1e-12"
      LABELS smoke)
endforeach()

# the same run on more ranks and more threads. These are not yet compared
# against each other; that needs the reduction output and comes next.
seissol_add_minitest(
    NAME smoke-box48-mpi4
    MESH box48
    PARAMETERS wave.par
    FILES material.yaml receiver.dat
    STEPS 20
    TIMESTEP 1e-4
    RANKS 4
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS "smoke;mpi")

seissol_add_minitest(
    NAME smoke-box48-threads4
    MESH box48
    PARAMETERS wave.par
    FILES material.yaml receiver.dat
    STEPS 20
    TIMESTEP 1e-4
    THREADS 4
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS "smoke;threads")

# plasticity is a runtime switch, so it needs no separate build
seissol_add_minitest(
    NAME smoke-skew48-plasticity
    MESH skew48
    PARAMETERS wave.par
    FILES material.yaml receiver.dat
    STEPS 20
    TIMESTEP 1e-4
    REQUIRE equations=elastic
    SET "equations.Plasticity=1"
    ACTIVITY "v1:1e-12"
    LABELS "smoke;plasticity")

seissol_add_minitest(
    NAME smoke-fault48
    MESH fault48
    PARAMETERS wave.par
    FILES material.yaml receiver.dat
    STEPS 20
    TIMESTEP 1e-4
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    LABELS "smoke;dr")


# ---------------------------------------------------------------------------
# seissol_add_comparison(NAME <name> REFERENCE <case> CANDIDATE <case>
#                        [TOLERANCE r] [PREFIX p] [LABELS label ...])
#
# Holds two runs against each other. A tolerance of zero demands identical
# output. The runner propagates a skip: if either input was skipped because
# this build lacks the configuration, so is the comparison.
# ---------------------------------------------------------------------------

function(seissol_add_comparison)
  cmake_parse_arguments(CMP
      ""
      "NAME;REFERENCE;CANDIDATE;TOLERANCE;PREFIX"
      "LABELS"
      ${ARGN})
  if (NOT CMP_TOLERANCE)
    set(CMP_TOLERANCE 0)
  endif()
  if (NOT CMP_PREFIX)
    set(CMP_PREFIX output/mini)
  endif()

  add_test(NAME verification.${CMP_NAME}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py compare
                   --name ${CMP_NAME}
                   --reference ${CMAKE_CURRENT_BINARY_DIR}/work/${CMP_REFERENCE}
                   --candidate ${CMAKE_CURRENT_BINARY_DIR}/work/${CMP_CANDIDATE}
                   --reference-result ${CMAKE_CURRENT_BINARY_DIR}/results/${CMP_REFERENCE}.json
                   --candidate-result ${CMAKE_CURRENT_BINARY_DIR}/results/${CMP_CANDIDATE}.json
                   --result ${CMAKE_CURRENT_BINARY_DIR}/results/${CMP_NAME}.json
                   --prefix ${CMP_PREFIX}
                   --tolerance ${CMP_TOLERANCE})
  set_tests_properties(verification.${CMP_NAME} PROPERTIES
      SKIP_RETURN_CODE 77
      FIXTURES_REQUIRED "case.${CMP_REFERENCE};case.${CMP_CANDIDATE};verification"
      TIMEOUT ${SEISSOL_MINITEST_TIMEOUT}
      LABELS "verification;${CMP_LABELS}")

  set_property(GLOBAL APPEND_STRING PROPERTY SEISSOL_DECLARED_CASES
      "    {\"name\": \"${CMP_NAME}\", \"mesh\": \"-\", \"requires\": []},\n")
endfunction()

# ---------------------------------------------------------------------------
# Stage 1: hold runs against each other. The receiver output carries about
# sixteen digits, so it can support a claim of identity; the analysis CSV is
# printed at the default stream precision and cannot.
# ---------------------------------------------------------------------------

seissol_add_comparison(
    NAME bitwise-box48-mpi4
    REFERENCE smoke-box48
    CANDIDATE smoke-box48-mpi4
    TOLERANCE 0
    LABELS "bitwise;mpi")

seissol_add_comparison(
    NAME bitwise-box48-threads4
    REFERENCE smoke-box48
    CANDIDATE smoke-box48-threads4
    TOLERANCE 0
    LABELS "bitwise;threads")

# ---------------------------------------------------------------------------
# Stage 2: the analytical error. A planar wave only satisfies a fully periodic
# mesh, so this is the mesh the measurement runs on. The thresholds are
# deliberately unset; see cases/thresholds.json.
# ---------------------------------------------------------------------------

seissol_add_minitest(
    NAME convergence-periodic162
    MESH periodic162
    PARAMETERS wave.par
    FILES material.yaml receiver.dat
    STEPS 100
    TIMESTEP 1e-4
    THRESHOLDS thresholds.json
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS convergence)

seissol_add_minitest(
    NAME convergence-periodicx72
    MESH periodicx72
    PARAMETERS wave.par
    FILES material.yaml receiver.dat
    STEPS 100
    TIMESTEP 1e-4
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS convergence)

# ---------------------------------------------------------------------------
# Coverage. Runs after the cases because it cleans up the fixture they
# require, and fails when a declared case has never been run by any build.
# ---------------------------------------------------------------------------

get_property(declared GLOBAL PROPERTY SEISSOL_DECLARED_CASES)
string(REGEX REPLACE ",\n$" "\n" declared "${declared}")
file(GENERATE
     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/declared.json
     CONTENT "{\n  \"cases\": [\n${declared}  ]\n}\n")

add_test(NAME verification.coverage
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py coverage
                 --declared ${CMAKE_CURRENT_BINARY_DIR}/declared.json
                 --results ${CMAKE_CURRENT_BINARY_DIR}/results)
set_tests_properties(verification.coverage PROPERTIES
                     FIXTURES_CLEANUP verification
                     LABELS "verification;coverage")
