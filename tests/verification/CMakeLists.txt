# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# ---------------------------------------------------------------------------
# Verification tests on generated mini meshes.
#
# The case list of Cases.cmake is the declared matrix. Every build registers
# all of it, for each of its configurations (or those SEISSOL_VERIFICATION_CONFIGS
# names), and skips what a configuration cannot run, so that a gap shows up as
# a skip instead of disappearing. The coverage test of a configuration turns
# those skips into a statement about what has been covered overall.
# ---------------------------------------------------------------------------

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# The mini meshes are written with numpy and h5py. Without them every mesh fails to generate, and
# every case on it is reported as not run, which says nothing about why; so it is said here.
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, h5py"
                RESULT_VARIABLE SEISSOL_MINIMESH_MODULES
                OUTPUT_QUIET ERROR_QUIET)
if (NOT SEISSOL_MINIMESH_MODULES EQUAL 0)
  message(WARNING "The verification tests need numpy and h5py for ${Python3_EXECUTABLE}. "
                  "Without them, no mini mesh can be generated and none of the cases runs.")
endif()

set(SEISSOL_MINITEST_TIMEOUT 300 CACHE STRING
    "Timeout in seconds for a single verification case")
set(SEISSOL_MINITEST_MPIEXEC_FLAGS "" CACHE STRING
    "Extra launcher flags for the verification cases, e.g. --oversubscribe on a \
machine with fewer cores than the rank counts the cases ask for")
# A case runs in a configuration of the executable, which the runner chooses with
# SEISSOL_CONFIGURATION; all of them take long in an unoptimized build (e.g. for the coverage).
set(SEISSOL_VERIFICATION_CONFIGS "all" CACHE STRING
    "Configurations to register the verification cases for: all of the build, first, or a list \
of their names")
set(SEISSOL_SNAPSHOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/references" CACHE PATH
    "Where the snapshot references live; may point to a separate repository")
# References are kept for few configurations, so that a deliberate change of
# the numerics can be accepted from a single build. Elastic carries by far the
# most cases, since rupture and plasticity need it. A build outside this list
# reports its snapshots as skipped, which the coverage test lists, so the gap
# stays visible; a case marked SNAPSHOT_EVERYWHERE is kept in every build.
set(SEISSOL_SNAPSHOT_CONFIGURATIONS "elastic-m0-o6-double-s1" CACHE STRING
    "Build configurations that keep snapshot references")
# The tolerances are a property of the environment a build runs in, relative
# to the one the references were computed in, not of the cases: identical
# environments reproduce the references bit for bit. Zero is therefore the
# default, and another environment sets its own tolerance from the deviation
# the failing comparison reports.
set(SEISSOL_SNAPSHOT_TOLERANCE_DOUBLE 0 CACHE STRING
    "Largest relative deviation from the snapshot references, double precision")
set(SEISSOL_SNAPSHOT_TOLERANCE_SINGLE 0 CACHE STRING
    "Largest relative deviation from the snapshot references, single precision")

# 2^-14 s. It lies below the CFL limit of every mini mesh up to order 8 for
# the materials in cases/, with a margin of at least 1.26, so the cap in
# FixTimeStep is what sets the step. Being a power of two, it also keeps
# N * dt and the accumulated time exact; the runner insists on both.
set(SEISSOL_MINITEST_TIMESTEP 6.103515625e-05)

# ---------------------------------------------------------------------------
# What this build can do.
#
# Only these axes are fixed at compile time. Plasticity, the friction law, the
# boundary conditions, the flux and the time stepping are all runtime
# parameters, so one build covers most of the matrix.
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

# not what the build can do, but where it runs; recorded with every snapshot
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Environment.json.in
    ${CMAKE_CURRENT_BINARY_DIR}/environment.json
    @ONLY
)

define_property(GLOBAL PROPERTY SEISSOL_DECLARED_CASES
                BRIEF_DOCS "declared verification cases"
                FULL_DOCS "declared verification cases")
define_property(GLOBAL PROPERTY SEISSOL_MINIMESHES
                BRIEF_DOCS "generated mini meshes"
                FULL_DOCS "generated mini meshes")

set(SEISSOL_MINIMESH_DIR ${CMAKE_CURRENT_BINARY_DIR}/meshes)

# ---------------------------------------------------------------------------
# seissol_add_minimesh(NAME <name> PRESET <preset> [ARGS <generator args>...])
#
# The generator's presets are topologies; which boundary condition sits on
# which plane is chosen here, per test, so the boundary condition matrix does
# not need a preset for every combination.
# ---------------------------------------------------------------------------

function(seissol_add_minimesh)
  cmake_parse_arguments(MM "" "NAME;PRESET" "ARGS" ${ARGN})
  add_test(NAME verification.mesh.${MM_NAME}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/minimesh.py generate ${MM_PRESET}
                   -o ${SEISSOL_MINIMESH_DIR}/${MM_NAME}.puml.h5
                   --manifest ${SEISSOL_MINIMESH_DIR}/${MM_NAME}.json
                   ${MM_ARGS})
  set_tests_properties(verification.mesh.${MM_NAME} PROPERTIES
                       FIXTURES_SETUP minimesh.${MM_NAME}
                       LABELS "verification;mesh")
  set_property(GLOBAL APPEND PROPERTY SEISSOL_MINIMESHES ${MM_NAME})
endfunction()

foreach (preset box6 box48 skew48 fault48 faultskew48 mat48 periodicx72 periodic162 open96)
  seissol_add_minimesh(NAME ${preset} PRESET ${preset})
endforeach()

# Dirichlet and Analytical are left out for now. Free surface with gravity is
# meant for acoustic media only. fsg2 is there for cells with two such faces,
# a case still to be looked at in the implementation. A Kuhn tetrahedron
# touches exactly one lower and one upper cube plane, on different axes, so
# two upper planes never share a cell; zmax with xmin gives one such cell.
seissol_add_minimesh(NAME box6-outflow PRESET box6 ARGS --bc all=outflow)
seissol_add_minimesh(NAME box6-mixed PRESET box6
                     ARGS --bc xmin=free-surface,xmax=outflow,ymin=outflow,ymax=free-surface,zmin=free-surface,zmax=outflow)
seissol_add_minimesh(NAME box6-fsg PRESET box6 ARGS --bc zmax=free-surface-gravity)
seissol_add_minimesh(NAME box6-fsg2 PRESET box6
                     ARGS --bc zmax=free-surface-gravity,xmin=free-surface-gravity)

add_test(NAME verification.minimesh.self-test
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/minimesh.py self-test)
set_tests_properties(verification.minimesh.self-test PROPERTIES
                     LABELS "verification;smoke")

# every receiver strictly inside a cell of every mesh; see minimesh.py locate
get_property(meshes GLOBAL PROPERTY SEISSOL_MINIMESHES)
set(mesh_files "")
set(mesh_fixtures "")
foreach (mesh ${meshes})
  list(APPEND mesh_files ${SEISSOL_MINIMESH_DIR}/${mesh}.puml.h5)
  list(APPEND mesh_fixtures minimesh.${mesh})
endforeach()
add_test(NAME verification.receivers
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/minimesh.py locate
                 ${CMAKE_CURRENT_SOURCE_DIR}/cases/receiver.dat ${mesh_files})
set_tests_properties(verification.receivers PROPERTIES
                     FIXTURES_REQUIRED "${mesh_fixtures}"
                     LABELS "verification;smoke")

add_test(NAME verification.faultreceivers
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/minimesh.py locate --fault
                 ${CMAKE_CURRENT_SOURCE_DIR}/cases/faultreceiver.dat
                 ${SEISSOL_MINIMESH_DIR}/fault48.puml.h5
                 ${SEISSOL_MINIMESH_DIR}/faultskew48.puml.h5)
set_tests_properties(verification.faultreceivers PROPERTIES
                     FIXTURES_REQUIRED "minimesh.fault48;minimesh.faultskew48"
                     LABELS "verification;smoke;dr")

# ---------------------------------------------------------------------------
# seissol_add_minitest(NAME <name> MESH <mesh> PARAMETERS <file>
#                      [MATERIAL file] [RANKS n] [THREADS n] [STEPS n]
#                      [TIMESTEP dt] [THRESHOLDS file] [EXPECT_ZERO]
#                      [EXPECT_FAILURE regex] [NO_SNAPSHOT] [SNAPSHOT_KEY key]
#                      [SNAPSHOT_EVERYWHERE]
#                      [REQUIRE key=value ...] [SET Section.Key=Value ...]
#                      [FILES name ...] [ACTIVITY column:threshold ...]
#                      [LABELS label ...])
#
# REQUIRE is matched against capabilities.json by the runner, not here: a case
# that this build cannot run has to remain a visible skip. File names and
# overrides may contain {equations}, {order} or {precision}; MATERIAL defaults
# to material-{equations}.yaml, and a build whose equation system has no such
# file yet skips the case, which keeps that gap visible as well.
#
# Every case that computes something gets a snapshot test holding it against
# a stored reference, one per numerical configuration as SNAPSHOT_KEY names
# it. Cases on several ranks or threads get none: their bit identity with the
# serial case is tested directly, and the serial case carries the snapshot.
# ---------------------------------------------------------------------------

function(seissol_add_minitest)
  cmake_parse_arguments(MT
      "EXPECT_ZERO;NO_SNAPSHOT;SNAPSHOT_EVERYWHERE"
      "NAME;MESH;PARAMETERS;MATERIAL;RANKS;THREADS;STEPS;TIMESTEP;THRESHOLDS;EXPECT_FAILURE;SNAPSHOT_KEY"
      "REQUIRE;SET;FILES;ACTIVITY;LABELS"
      ${ARGN})
  if (NOT MT_NAME OR NOT MT_MESH OR NOT MT_PARAMETERS)
    message(FATAL_ERROR "seissol_add_minitest needs NAME, MESH and PARAMETERS")
  endif()
  get_property(meshes GLOBAL PROPERTY SEISSOL_MINIMESHES)
  if (NOT MT_MESH IN_LIST meshes)
    message(FATAL_ERROR "unknown mini mesh '${MT_MESH}' in case '${MT_NAME}'")
  endif()
  if (NOT MT_RANKS)
    set(MT_RANKS 1)
  endif()
  if (NOT MT_THREADS)
    set(MT_THREADS 1)
  endif()
  if (NOT MT_TIMESTEP)
    set(MT_TIMESTEP ${SEISSOL_MINITEST_TIMESTEP})
  endif()
  if (NOT MT_MATERIAL)
    set(MT_MATERIAL "material-{equations}.yaml")
  endif()
  list(APPEND MT_FILES ${MT_MATERIAL} receiver.dat)
  list(APPEND MT_SET "equations.MaterialFileName='${MT_MATERIAL}'")

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
      --seissol-configuration ${SEISSOL_VERIFICATION_CONFIG}
      --capabilities ${SEISSOL_VERIFICATION_DIR}/capabilities.json
      --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/cases
      --work-dir ${SEISSOL_VERIFICATION_DIR}/work/${MT_NAME}
      --result ${SEISSOL_VERIFICATION_DIR}/results/${MT_NAME}.json
      --mesh ${SEISSOL_MINIMESH_DIR}/${MT_MESH}.puml.h5
      --mesh-manifest ${SEISSOL_MINIMESH_DIR}/${MT_MESH}.json
      --parameters ${MT_PARAMETERS}
      --ranks ${MT_RANKS}
      --threads ${MT_THREADS}
      --timestep ${MT_TIMESTEP})
  if (MT_STEPS)
    list(APPEND arguments --steps ${MT_STEPS})
  endif()
  if (MT_THRESHOLDS)
    list(APPEND arguments --thresholds ${CMAKE_CURRENT_SOURCE_DIR}/cases/${MT_THRESHOLDS})
  endif()
  if (MT_EXPECT_ZERO)
    list(APPEND arguments --expect-zero)
  endif()
  if (MT_EXPECT_FAILURE)
    list(APPEND arguments "--expect-failure=${MT_EXPECT_FAILURE}")
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
  elseif (TESTING_COMMAND)
    # a single rank runs under TESTING_COMMAND, as the unit and smoke tests do: a runtime that
    # needs a launcher even then, like the one of NVHPC, is given one there
    separate_arguments(launcher NATIVE_COMMAND "${TESTING_COMMAND}")
    foreach (token ${launcher})
      list(APPEND arguments "--launcher=${token}")
    endforeach()
  endif()

  add_test(NAME verification.${MT_NAME}${SEISSOL_VERIFICATION_SUFFIX}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py ${arguments})

  math(EXPR slots "${MT_RANKS} * ${MT_THREADS}")
  set_tests_properties(verification.${MT_NAME}${SEISSOL_VERIFICATION_SUFFIX} PROPERTIES
      # the runner reports a configuration this build lacks as skipped, so that
      # it is neither a failure nor silently absent
      SKIP_RETURN_CODE 77
      FIXTURES_SETUP case.${SEISSOL_VERIFICATION_CONFIG}.${MT_NAME}
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

  if (NOT MT_SNAPSHOT_KEY)
    set(MT_SNAPSHOT_KEY "{equations}-m{mechanisms}-o{order}-{precision}-s{fused}")
  endif()
  set(snapshot_scope "")
  if (NOT MT_SNAPSHOT_EVERYWHERE)
    set(snapshot_scope --only-configurations ${SEISSOL_SNAPSHOT_CONFIGURATIONS})
  endif()
  if (NOT MT_NO_SNAPSHOT AND NOT MT_EXPECT_ZERO AND NOT MT_EXPECT_FAILURE
      AND MT_RANKS EQUAL 1 AND MT_THREADS EQUAL 1)
    add_test(NAME verification.snapshot.${MT_NAME}${SEISSOL_VERIFICATION_SUFFIX}
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py snapshot
                     --name ${MT_NAME}
                     --case-result ${SEISSOL_VERIFICATION_DIR}/results/${MT_NAME}.json
                     --work ${SEISSOL_VERIFICATION_DIR}/work/${MT_NAME}
                     --reference-dir ${SEISSOL_SNAPSHOT_DIR}
                     --key ${MT_SNAPSHOT_KEY}
                     --capabilities ${SEISSOL_VERIFICATION_DIR}/capabilities.json
                     --environment ${CMAKE_CURRENT_BINARY_DIR}/environment.json
                     --source-dir ${PROJECT_SOURCE_DIR}
                     --configuration ${SEISSOL_BUILD_CONFIGURATION}
                     ${snapshot_scope}
                     --tolerance-double ${SEISSOL_SNAPSHOT_TOLERANCE_DOUBLE}
                     --tolerance-single ${SEISSOL_SNAPSHOT_TOLERANCE_SINGLE}
                     --result ${SEISSOL_VERIFICATION_DIR}/results/snapshot-${MT_NAME}.json)
    # with the labels of its case as well, so that a selection by label takes or leaves both
    set(snapshot_labels verification snapshot ${MT_LABELS})
    set_tests_properties(verification.snapshot.${MT_NAME}${SEISSOL_VERIFICATION_SUFFIX} PROPERTIES
        SKIP_RETURN_CODE 77
        FIXTURES_REQUIRED "case.${SEISSOL_VERIFICATION_CONFIG}.${MT_NAME};verification"
        TIMEOUT ${SEISSOL_MINITEST_TIMEOUT}
        LABELS "${snapshot_labels}")
    set_property(GLOBAL APPEND_STRING PROPERTY SEISSOL_DECLARED_CASES
        "    {\"name\": \"snapshot-${MT_NAME}\", \"mesh\": \"${MT_MESH}\", \"requires\": []},\n")
  endif()
endfunction()

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

  set(directory ${SEISSOL_VERIFICATION_DIR})
  add_test(NAME verification.${CMP_NAME}${SEISSOL_VERIFICATION_SUFFIX}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py compare
                   --name ${CMP_NAME}
                   --reference ${directory}/work/${CMP_REFERENCE}
                   --candidate ${directory}/work/${CMP_CANDIDATE}
                   --reference-result ${directory}/results/${CMP_REFERENCE}.json
                   --candidate-result ${directory}/results/${CMP_CANDIDATE}.json
                   --result ${directory}/results/${CMP_NAME}.json
                   --prefix ${CMP_PREFIX}
                   --tolerance ${CMP_TOLERANCE})
  set(fixtures
      case.${SEISSOL_VERIFICATION_CONFIG}.${CMP_REFERENCE}
      case.${SEISSOL_VERIFICATION_CONFIG}.${CMP_CANDIDATE}
      verification)
  set_tests_properties(verification.${CMP_NAME}${SEISSOL_VERIFICATION_SUFFIX} PROPERTIES
      SKIP_RETURN_CODE 77
      FIXTURES_REQUIRED "${fixtures}"
      TIMEOUT ${SEISSOL_MINITEST_TIMEOUT}
      LABELS "verification;${CMP_LABELS}")

  set_property(GLOBAL APPEND_STRING PROPERTY SEISSOL_DECLARED_CASES
      "    {\"name\": \"${CMP_NAME}\", \"mesh\": \"-\", \"requires\": []},\n")
endfunction()

# ---------------------------------------------------------------------------
# seissol_add_rupture(NAME <name> FL <law> FAULT <file> [SET ...] [MESH mesh]
#                     [RANKS n] [EXPECT_FAILURE regex] [EXPECT_ZERO])
#
# A case of a friction law on fault48 (see the friction laws in Cases.cmake).
# ---------------------------------------------------------------------------

function(seissol_add_rupture)
  cmake_parse_arguments(RUP
      "EXPECT_ZERO;NO_SNAPSHOT;SNAPSHOT_EVERYWHERE"
      "NAME;FL;FAULT;MESH;RANKS;EXPECT_FAILURE"
      "SET;LABELS"
      ${ARGN})
  if (NOT RUP_MESH)
    set(RUP_MESH fault48)
  endif()
  set(options "")
  if (RUP_EXPECT_ZERO)
    list(APPEND options EXPECT_ZERO)
  else()
    list(APPEND options ACTIVITY "SRs|SRd:1e-6")
  endif()
  if (RUP_EXPECT_FAILURE)
    list(APPEND options EXPECT_FAILURE "${RUP_EXPECT_FAILURE}")
  endif()
  if (RUP_RANKS)
    list(APPEND options RANKS ${RUP_RANKS})
  endif()
  if (RUP_NO_SNAPSHOT)
    list(APPEND options NO_SNAPSHOT)
  endif()
  seissol_add_minitest(
      NAME ${RUP_NAME}
      MESH ${RUP_MESH}
      PARAMETERS rupture.par
      # the quadrature rule on the fault changes the result of every case with one
      SNAPSHOT_KEY "{equations}-m{mechanisms}-o{order}-{precision}-s{fused}-{dr-quad-rule}"
      FILES ${RUP_FAULT} faultreceiver.dat
      STEPS 64
      REQUIRE equations=elastic
      SET "DynamicRupture.FL=${RUP_FL}"
          "DynamicRupture.ModelFileName='${RUP_FAULT}'"
          ${RUP_SET}
      ${options}
      LABELS "dr;${RUP_LABELS}")
endfunction()

# ---------------------------------------------------------------------------
# The cases of Cases.cmake, for every configuration they are registered for.
#
# The tests of a configuration carry its name, as in
# verification.smoke-box48<elastic-linearck-o6-f64-stroud>, and work in a
# directory of their own; the runner chooses the configuration of the
# executable with SEISSOL_CONFIGURATION. The coverage test of a configuration
# runs after its cases because it cleans up the fixture they require, and fails
# when a declared case has never been run by any build.
# ---------------------------------------------------------------------------

if (SEISSOL_VERIFICATION_CONFIGS STREQUAL "all")
  set(verification_configs ${SEISSOL_CONFIG_NAMES})
elseif (SEISSOL_VERIFICATION_CONFIGS STREQUAL "first")
  list(GET SEISSOL_CONFIG_NAMES 0 verification_configs)
else()
  set(verification_configs ${SEISSOL_VERIFICATION_CONFIGS})
endif()
message(STATUS "Verification: cases for ${verification_configs}")

foreach (SEISSOL_VERIFICATION_CONFIG IN LISTS verification_configs)
  list(FIND SEISSOL_CONFIG_NAMES ${SEISSOL_VERIFICATION_CONFIG} index)
  if (index LESS 0)
    message(FATAL_ERROR "SEISSOL_VERIFICATION_CONFIGS: ${SEISSOL_VERIFICATION_CONFIG} is not "
                        "built; the build holds ${SEISSOL_CONFIG_NAMES}")
  endif()
  # what the configuration is, for capabilities.json and the snapshot keys
  foreach (property MATERIALS SOLVERS MECHANISMS ORDERS PRECISIONS DRQUADRULES SIMULATIONS)
    list(GET SEISSOL_CONFIG_${property} ${index} SEISSOL_CAP_${property})
  endforeach()
  set(SEISSOL_BUILD_CONFIGURATION
      "${SEISSOL_CAP_MATERIALS}-m${SEISSOL_CAP_MECHANISMS}-o${SEISSOL_CAP_ORDERS}-${SEISSOL_CAP_PRECISIONS}-s${SEISSOL_CAP_SIMULATIONS}")
  set(SEISSOL_VERIFICATION_SUFFIX "<${SEISSOL_VERIFICATION_CONFIG}>")
  set(SEISSOL_VERIFICATION_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SEISSOL_VERIFICATION_CONFIG})
  configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/Capabilities.json.in
      ${SEISSOL_VERIFICATION_DIR}/capabilities.json
      @ONLY
  )

  set_property(GLOBAL PROPERTY SEISSOL_DECLARED_CASES "")
  include(${CMAKE_CURRENT_SOURCE_DIR}/Cases.cmake)

  get_property(declared GLOBAL PROPERTY SEISSOL_DECLARED_CASES)
  string(REGEX REPLACE ",\n$" "\n" declared "${declared}")
  file(GENERATE
       OUTPUT ${SEISSOL_VERIFICATION_DIR}/declared.json
       CONTENT "{\n  \"cases\": [\n${declared}  ]\n}\n")

  add_test(NAME verification.coverage${SEISSOL_VERIFICATION_SUFFIX}
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py coverage
                   --declared ${SEISSOL_VERIFICATION_DIR}/declared.json
                   --results ${SEISSOL_VERIFICATION_DIR}/results)
  set_tests_properties(verification.coverage${SEISSOL_VERIFICATION_SUFFIX} PROPERTIES
                       FIXTURES_CLEANUP verification
                       LABELS "verification;coverage")
endforeach()
