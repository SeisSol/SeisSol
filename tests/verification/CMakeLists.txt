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

set(SEISSOL_BUILD_CONFIGURATION
    "${EQUATIONS}-m${NUMBER_OF_MECHANISMS}-o${ORDER}-${PRECISION}-s${NUMBER_OF_FUSED_SIMULATIONS}")
message(STATUS "Verification: build configuration ${SEISSOL_BUILD_CONFIGURATION}")

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Capabilities.json.in
    ${CMAKE_CURRENT_BINARY_DIR}/capabilities.json
    @ONLY
)
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
      --capabilities ${CMAKE_CURRENT_BINARY_DIR}/capabilities.json
      --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/cases
      --work-dir ${CMAKE_CURRENT_BINARY_DIR}/work/${MT_NAME}
      --result ${CMAKE_CURRENT_BINARY_DIR}/results/${MT_NAME}.json
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

  if (NOT MT_SNAPSHOT_KEY)
    set(MT_SNAPSHOT_KEY "{equations}-m{mechanisms}-o{order}-{precision}-s{fused}")
  endif()
  set(snapshot_scope "")
  if (NOT MT_SNAPSHOT_EVERYWHERE)
    set(snapshot_scope --only-configurations ${SEISSOL_SNAPSHOT_CONFIGURATIONS})
  endif()
  if (NOT MT_NO_SNAPSHOT AND NOT MT_EXPECT_ZERO AND NOT MT_EXPECT_FAILURE
      AND MT_RANKS EQUAL 1 AND MT_THREADS EQUAL 1)
    add_test(NAME verification.snapshot.${MT_NAME}
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py snapshot
                     --name ${MT_NAME}
                     --case-result ${CMAKE_CURRENT_BINARY_DIR}/results/${MT_NAME}.json
                     --work ${CMAKE_CURRENT_BINARY_DIR}/work/${MT_NAME}
                     --reference-dir ${SEISSOL_SNAPSHOT_DIR}
                     --key ${MT_SNAPSHOT_KEY}
                     --capabilities ${CMAKE_CURRENT_BINARY_DIR}/capabilities.json
                     --environment ${CMAKE_CURRENT_BINARY_DIR}/environment.json
                     --source-dir ${PROJECT_SOURCE_DIR}
                     --configuration ${SEISSOL_BUILD_CONFIGURATION}
                     ${snapshot_scope}
                     --tolerance-double ${SEISSOL_SNAPSHOT_TOLERANCE_DOUBLE}
                     --tolerance-single ${SEISSOL_SNAPSHOT_TOLERANCE_SINGLE}
                     --result ${CMAKE_CURRENT_BINARY_DIR}/results/snapshot-${MT_NAME}.json)
    set_tests_properties(verification.snapshot.${MT_NAME} PROPERTIES
        SKIP_RETURN_CODE 77
        FIXTURES_REQUIRED "case.${MT_NAME};verification"
        TIMEOUT ${SEISSOL_MINITEST_TIMEOUT}
        LABELS "verification;snapshot")
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
# Stage 0: does it run at all, and did it compute anything.
# ---------------------------------------------------------------------------

set(SEISSOL_WAVE_MATERIALS "equations=elastic|anisotropic|viscoelastic2|poroelastic|acoustic")

foreach (mesh box48 skew48)
  seissol_add_minitest(
      NAME smoke-${mesh}
      MESH ${mesh}
      PARAMETERS wave.par
      STEPS 20
      REQUIRE ${SEISSOL_WAVE_MATERIALS}
      ACTIVITY "v1:1e-12"
      LABELS smoke)
endforeach()

seissol_add_minitest(
    NAME smoke-box48-mpi4
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    RANKS 4
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS "smoke;mpi")

seissol_add_minitest(
    NAME smoke-box48-threads4
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    THREADS 4
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS "smoke;threads")

# plasticity is a runtime switch, so it needs no separate build; its
# parameters come from the material file
seissol_add_minitest(
    NAME smoke-skew48-plasticity
    MESH skew48
    PARAMETERS wave.par
    MATERIAL material-plastic.yaml
    STEPS 20
    REQUIRE equations=elastic
    SET "equations.Plasticity=1"
    ACTIVITY "v1:1e-12"
    LABELS "smoke;plasticity")

# without a friction law the fault is welded: NoFault passes the traction of
# the continuous-contact state through unchanged
seissol_add_minitest(
    NAME smoke-fault48
    MESH fault48
    PARAMETERS wave.par
    STEPS 20
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    LABELS "smoke;dr")

# ---------------------------------------------------------------------------
# Stage 1: invariants that need no reference data.
#
# A zero state without sources has to stay exactly zero. Every linear operator
# maps zero to zero whatever its sign or rotation, so this does not test the
# mathematics of a boundary condition; it tests that nothing is injected and
# that every operator is finite, since 0 * NaN is NaN. That covers
# uninitialised buffers and halos, material parameters missing in some cell,
# boundary data that is not zero, and nonlinear updates that act on a state
# they should leave alone.
# ---------------------------------------------------------------------------

set(SEISSOL_ZERO_STATE "IniCondition.cICType='Zero'")

foreach (mesh box6 box6-outflow box6-mixed box48 skew48 periodic162 mat48)
  seissol_add_minitest(
      NAME null-${mesh}
      MESH ${mesh}
      PARAMETERS wave.par
      STEPS 20
      SET ${SEISSOL_ZERO_STATE}
      REQUIRE ${SEISSOL_WAVE_MATERIALS}
      EXPECT_ZERO
      LABELS null)
endforeach()

# halos and the threaded paths are where uninitialised memory hides
seissol_add_minitest(
    NAME null-skew48-mpi4
    MESH skew48
    PARAMETERS wave.par
    STEPS 20
    RANKS 4
    SET ${SEISSOL_ZERO_STATE}
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    EXPECT_ZERO
    LABELS "null;mpi")

seissol_add_minitest(
    NAME null-skew48-threads4
    MESH skew48
    PARAMETERS wave.par
    STEPS 20
    THREADS 4
    SET ${SEISSOL_ZERO_STATE}
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    EXPECT_ZERO
    LABELS "null;threads")

foreach (mesh box6-fsg box6-fsg2)
  seissol_add_minitest(
      NAME null-${mesh}
      MESH ${mesh}
      PARAMETERS wave.par
      STEPS 20
      SET ${SEISSOL_ZERO_STATE}
      REQUIRE equations=acoustic
      EXPECT_ZERO
      LABELS "null;fsg")
endforeach()

# the initial stress is hydrostatic and below yield, so the plasticity update
# must leave a zero perturbation untouched
seissol_add_minitest(
    NAME null-skew48-plasticity
    MESH skew48
    PARAMETERS wave.par
    MATERIAL material-plastic.yaml
    STEPS 20
    SET ${SEISSOL_ZERO_STATE} "equations.Plasticity=1"
    REQUIRE equations=elastic
    EXPECT_ZERO
    LABELS "null;plasticity")

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

# A welded fault is the absence of a fault: fault48 and box48 are the same
# mesh apart from the tag on the plane x = 0, and the rupture quadrature is
# exact for the polynomial degrees involved, so the two runs differ by
# rounding only. A sign or rotation error in the rupture flux gives a
# difference of order one; the bound sits between the two regimes and gets
# tightened once the actual difference has been measured.
seissol_add_comparison(
    NAME identity-fault48-nofault
    REFERENCE smoke-box48
    CANDIDATE smoke-fault48
    TOLERANCE 1e-6
    LABELS "identity;dr")


# ---------------------------------------------------------------------------
# Friction laws.
#
# seissol_add_rupture(NAME <name> FL <law> FAULT <file> [SET ...] [MESH mesh]
#                     [RANKS n] [EXPECT_FAILURE regex] [EXPECT_ZERO])
#
# Every law runs on fault48 from rest, with a nucleation stress that breaks
# the fault within the run. The activity demand is a slip rate above 1e-6 m/s
# in either fault direction: ten orders of magnitude above the creep rate the
# rate-and-state laws start from, and exactly zero for a locked fault.
#
# Which easi file a law gets matters beyond its own parameters: SeisSol
# chooses between the Cartesian and the traction parametrization, and
# switches the forced rupture time on, by what the file provides. Each law
# therefore gets a file with exactly its parameters.
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
      # the rupture method changes the result of every case with a fault
      SNAPSHOT_KEY "{equations}-m{mechanisms}-o{order}-{precision}-s{fused}-{rupture-method}"
      FILES ${RUP_FAULT} faultreceiver.dat
      STEPS 64
      REQUIRE equations=elastic
      SET "DynamicRupture.FL=${RUP_FL}"
          "DynamicRupture.ModelFileName='${RUP_FAULT}'"
          ${RUP_SET}
      ${options}
      LABELS "dr;${RUP_LABELS}")
endfunction()

set(SEISSOL_RS
    "DynamicRupture.RS_f0=0.6"
    "DynamicRupture.RS_sr0=1e-6"
    "DynamicRupture.RS_b=0.012"
    "DynamicRupture.RS_iniSlipRate1=1e-16"
    "DynamicRupture.RS_iniSlipRate2=0.0")
set(SEISSOL_RS_FVW
    "DynamicRupture.RS_f0=0.6"
    "DynamicRupture.RS_sr0=1e-6"
    "DynamicRupture.RS_b=0.014"
    "DynamicRupture.RS_muW=0.2"
    "DynamicRupture.RS_iniSlipRate1=1e-16"
    "DynamicRupture.RS_iniSlipRate2=0.0")
set(SEISSOL_TP
    "DynamicRupture.thermalPress=1"
    "DynamicRupture.tp_iniTemp=483.15"
    "DynamicRupture.tp_iniPressure=0.0"
    "DynamicRupture.tp_thermalDiffusivity=1.0e-6"
    "DynamicRupture.tp_heatCapacity=2.7e6"
    "DynamicRupture.tp_undrainedTPResponse=0.1e6")
# delays the nucleation stress past the end of the run
set(SEISSOL_NO_NUCLEATION "DynamicRupture.s_0=1.0")

seissol_add_rupture(NAME rupture-lsw FL 16 FAULT fault-lsw.yaml)
# the forced rupture time is a separate branch, taken only if the file
# provides it; without nucleation it is the only thing breaking the fault
seissol_add_rupture(NAME rupture-lsw-forced FL 16 FAULT fault-lsw-forced.yaml
                    SET ${SEISSOL_NO_NUCLEATION})
seissol_add_rupture(NAME rupture-lsw-bimaterial FL 6 FAULT fault-lsw.yaml
                    SET "DynamicRupture.pc_vstar=0.1" "DynamicRupture.pc_prakashLength=0.1")
seissol_add_rupture(NAME rupture-lsw-tpapprox FL 1058 FAULT fault-lsw.yaml)
# FL=2 is documented as a deprecated spelling of FL=16; its identity with
# FL=16 is tested directly, so it needs no snapshot of its own
seissol_add_rupture(NAME rupture-lsw-legacy FL 2 FAULT fault-lsw.yaml NO_SNAPSHOT)

foreach (law aging:3 slip:4 svw:7)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  seissol_add_rupture(NAME rupture-rs-${name} FL ${number} FAULT fault-rs.yaml
                      SET ${SEISSOL_RS})
  seissol_add_rupture(NAME rupture-rs-${name}-tp FL ${number} FAULT fault-rs-tp.yaml
                      SET ${SEISSOL_RS} ${SEISSOL_TP})
endforeach()
seissol_add_rupture(NAME rupture-rs-fvw FL 103 FAULT fault-rs-fvw.yaml
                    SET ${SEISSOL_RS_FVW})
seissol_add_rupture(NAME rupture-rs-fvw-tp FL 103 FAULT fault-rs-fvw-tp.yaml
                    SET ${SEISSOL_RS_FVW} ${SEISSOL_TP})

foreach (law yoffe:33 gaussian:34 delta:35)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  seissol_add_rupture(NAME rupture-imposed-${name} FL ${number} FAULT fault-imposed.yaml
                      SET "DynamicRupture.SlipRateOutputType=0")
endforeach()

# disabled on purpose in SeisSol; it has to be refused at initialization
seissol_add_rupture(NAME rupture-rs-aging-nucleation FL 101 FAULT fault-rs.yaml
                    SET ${SEISSOL_RS}
                    EXPECT_FAILURE "friction law 101 currently disabled")

# a prestressed fault below its strength is locked and must not radiate: the
# initial stress is background, and nothing of it may leak into the wavefield
seissol_add_rupture(NAME null-faultskew48-lsw FL 16 FAULT fault-lsw.yaml
                    MESH faultskew48 SET ${SEISSOL_NO_NUCLEATION} EXPECT_ZERO
                    LABELS null)

# rupture across rank boundaries: fault faces end up in the copy layer
seissol_add_rupture(NAME rupture-lsw-mpi4 FL 16 FAULT fault-lsw.yaml RANKS 4
                    LABELS mpi)
seissol_add_comparison(
    NAME bitwise-rupture-lsw-mpi4
    REFERENCE rupture-lsw
    CANDIDATE rupture-lsw-mpi4
    TOLERANCE 0
    LABELS "bitwise;mpi;dr")

# a deprecated spelling of the same law has to give the same run
seissol_add_comparison(
    NAME identity-lsw-legacy
    REFERENCE rupture-lsw
    CANDIDATE rupture-lsw-legacy
    TOLERANCE 0
    LABELS "identity;dr")

# ---------------------------------------------------------------------------
# Stage 2: the analytical error. A planar wave only satisfies a fully periodic
# mesh, so this is the mesh the measurement runs on. 256 steps of 2^-14 s let
# the P wave cross about one wavelength. The thresholds are deliberately
# unset; see cases/thresholds.json.
# ---------------------------------------------------------------------------

seissol_add_minitest(
    NAME convergence-periodic162
    MESH periodic162
    PARAMETERS wave.par
    STEPS 256
    THRESHOLDS thresholds.json
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS convergence)

seissol_add_minitest(
    NAME convergence-periodicx72
    MESH periodicx72
    PARAMETERS wave.par
    STEPS 256
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
