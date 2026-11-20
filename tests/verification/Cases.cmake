# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The verification cases; CMakeLists.txt includes this file once for every configuration it
# registers them for, after setting what describes the configuration.

# ---------------------------------------------------------------------------
# Stage 0: does it run at all, and did it compute anything.
# ---------------------------------------------------------------------------

set(SEISSOL_WAVE_MATERIALS "equations=elastic|anisotropic|viscoelastic|poroelastic|acoustic")

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

# The task graph on the host with a single thread per rank, through two epochs (the receivers
# write after 10 of the 20 steps): each epoch is a parallel region of one thread, whose halo
# exchange runs in detached tasks. LLVM's OpenMP runtime fails at the second one before version 19;
# SeisSol keeps the static schedule there.
seissol_add_minitest(
    NAME smoke-box48-taskgraph-mpi2
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    RANKS 2
    SET "Output.ReceiverOutputInterval=0.0006103515625"
    ENV "SEISSOL_HOST_TASK_GRAPH=1"
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS "smoke;mpi;taskgraph")

# plasticity is a runtime switch, so it needs no separate build. Under the
# parameters of material-plastic.yaml, the wave of the case stays elastic: it
# reaches less than 1 Pa, where the material yields at about 5 MPa. The case
# takes those of plasticity-weak.yaml, under which it yields, and the plastic
# strain relaxes within a step, so that yielding shows.
seissol_add_minitest(
    NAME smoke-skew48-plasticity
    MESH skew48
    PARAMETERS wave.par
    MATERIAL material-plastic.yaml
    STEPS 20
    REQUIRE equations=elastic
    SET "equations.Plasticity=1" "equations.Tv=6.103515625e-05"
        "equations.PlasticityFileName='plasticity-weak.yaml'"
    FILES plasticity-weak.yaml
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

# the Dirichlet faces drive a wave into the body at rest; their setup takes the
# nodes of the faces, which the fused configurations lay out differently
seissol_add_minitest(
    NAME smoke-box48-dirichlet
    MESH box48-dirichlet
    PARAMETERS wave.par
    STEPS 20
    SET "IniCondition.cICType='Zero'" "equations.BoundaryFileName='boundary-v1.yaml'"
    FILES boundary-v1.yaml
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    LABELS smoke)

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
# sixteen digits, so it can support a claim of identity.
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

seissol_add_comparison(
    NAME bitwise-box48-taskgraph-mpi2
    REFERENCE smoke-box48
    CANDIDATE smoke-box48-taskgraph-mpi2
    TOLERANCE 0
    LABELS "bitwise;mpi;taskgraph")

# Along the time stepping plan (SEISSOL_TIMESTEPPING_PLAN, see TimeManager::followPlan), with the
# ghost clusters of two ranks, through two epochs: the same run, bit for bit.
seissol_add_minitest(
    NAME smoke-box48-plan-mpi2
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    RANKS 2
    SET "Output.ReceiverOutputInterval=0.0006103515625"
    ENV "SEISSOL_TIMESTEPPING_PLAN=1"
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS "smoke;mpi;plan")

seissol_add_comparison(
    NAME bitwise-box48-plan-mpi2
    REFERENCE smoke-box48
    CANDIDATE smoke-box48-plan-mpi2
    TOLERANCE 0
    LABELS "bitwise;mpi;plan")

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
# Stage 1: simulation copies. Twice as many simulations as the configuration
# fuses are held by two copies of every cell. Each simulation has to carry a
# wave of its own, a zero state has to stay zero in every copy, and the copies
# exchange their halos as the cells do.
# ---------------------------------------------------------------------------

math(EXPR SEISSOL_COPIES_SIMULATIONS "2 * ${SEISSOL_CAP_SIMULATIONS}")

seissol_add_minitest(
    NAME smoke-box48-copies
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    SET "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}"
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;copies")

seissol_add_minitest(
    NAME smoke-box48-copies-mpi4
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    RANKS 4
    SET "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}"
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    ACTIVITY "v1:1e-12"
    LABELS "smoke;mpi;copies")

seissol_add_minitest(
    NAME null-skew48-copies-mpi4
    MESH skew48
    PARAMETERS wave.par
    STEPS 20
    RANKS 4
    SET ${SEISSOL_ZERO_STATE} "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}"
    REQUIRE ${SEISSOL_WAVE_MATERIALS}
    EXPECT_ZERO
    LABELS "null;mpi;copies")

seissol_add_comparison(
    NAME bitwise-box48-copies-mpi4
    REFERENCE smoke-box48-copies
    CANDIDATE smoke-box48-copies-mpi4
    TOLERANCE 0
    LABELS "bitwise;mpi;copies")

# The copies hold their simulations as the fused simulations of a cell hold
# theirs. The waves of a run are offset in phase by 2 pi s / S for its S
# simulations, so the simulation 2s of the run of the copies carries the wave
# of the simulation s of the run of the fused ones only, whichever copy and
# which fused simulation of it holds it.
set(SEISSOL_COPIES_OF_FUSED "")
math(EXPR _copies_last "${SEISSOL_CAP_SIMULATIONS} - 1")
foreach (_copies_simulation RANGE ${_copies_last})
  math(EXPR _copies_candidate "2 * ${_copies_simulation}")
  list(APPEND SEISSOL_COPIES_OF_FUSED "${_copies_simulation}:${_copies_candidate}")
endforeach()
seissol_add_comparison(
    NAME bitwise-box48-copies-fused
    REFERENCE smoke-box48
    CANDIDATE smoke-box48-copies
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_OF_FUSED}
    LABELS "bitwise;copies")

# A copy of the cells has the material of the simulations it holds. In a run
# whose second copy takes another material, the simulations of each copy are
# those of the run whose copies all take its material: their waves are set up
# with it, travel at its speeds, and meet the cells of their own copy only.
seissol_add_minitest(
    NAME smoke-box48-copies-soft
    MESH box48
    PARAMETERS wave.par
    MATERIAL material-elastic-soft.yaml
    STEPS 20
    SET "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}"
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;copies")
set(SEISSOL_COPIES_MATERIALS "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}")
set(SEISSOL_COPIES_FIRST "")
set(SEISSOL_COPIES_SECOND "")
math(EXPR _copies_last "${SEISSOL_COPIES_SIMULATIONS} - 1")
foreach (_copies_simulation RANGE ${_copies_last})
  if (_copies_simulation LESS SEISSOL_CAP_SIMULATIONS)
    list(APPEND SEISSOL_COPIES_FIRST "${_copies_simulation}:${_copies_simulation}")
  else()
    list(APPEND SEISSOL_COPIES_MATERIALS
         "equations.MaterialFileName${_copies_simulation}='material-elastic-soft.yaml'")
    list(APPEND SEISSOL_COPIES_SECOND "${_copies_simulation}:${_copies_simulation}")
  endif()
endforeach()
seissol_add_minitest(
    NAME smoke-box48-copies-materials
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    SET ${SEISSOL_COPIES_MATERIALS}
    FILES material-elastic-soft.yaml
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;copies")
seissol_add_comparison(
    NAME bitwise-box48-copies-materials-first
    REFERENCE smoke-box48-copies
    CANDIDATE smoke-box48-copies-materials
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_FIRST}
    LABELS "bitwise;copies")
seissol_add_comparison(
    NAME bitwise-box48-copies-materials-second
    REFERENCE smoke-box48-copies-soft
    CANDIDATE smoke-box48-copies-materials
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_SECOND}
    LABELS "bitwise;copies")

# An easi initial condition per simulation: every second simulation starts
# from the other state, in both copies of the cells and in the fused
# simulations of each. The runs whose simulations all start from the one or
# the other state hold them alike.
set(SEISSOL_EASI_STATE "IniCondition.cICType='easi'" "IniCondition.hasTime=0"
    "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}")
foreach (_copies_axis x y)
  seissol_add_minitest(
      NAME easi-box48-copies-${_copies_axis}
      MESH box48
      PARAMETERS wave.par
      STEPS 20
      SET ${SEISSOL_EASI_STATE} "IniCondition.filename='inicondition-${_copies_axis}.yaml'"
      FILES inicondition-${_copies_axis}.yaml
      REQUIRE equations=elastic
      ACTIVITY "v1|v2:1e-12"
      NO_SNAPSHOT
      LABELS "smoke;copies")
endforeach()
set(SEISSOL_COPIES_STATES ${SEISSOL_EASI_STATE} "IniCondition.filename='inicondition-x.yaml'")
set(SEISSOL_COPIES_EVEN "")
set(SEISSOL_COPIES_ODD "")
math(EXPR _copies_last "${SEISSOL_COPIES_SIMULATIONS} - 1")
foreach (_copies_simulation RANGE ${_copies_last})
  math(EXPR _copies_odd "${_copies_simulation} % 2")
  if (_copies_odd)
    list(APPEND SEISSOL_COPIES_STATES
         "IniCondition.filename${_copies_simulation}='inicondition-y.yaml'")
    list(APPEND SEISSOL_COPIES_ODD "${_copies_simulation}:${_copies_simulation}")
  else()
    list(APPEND SEISSOL_COPIES_EVEN "${_copies_simulation}:${_copies_simulation}")
  endif()
endforeach()
seissol_add_minitest(
    NAME easi-box48-copies-states
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    SET ${SEISSOL_COPIES_STATES}
    FILES inicondition-x.yaml inicondition-y.yaml
    REQUIRE equations=elastic
    ACTIVITY "v1|v2:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;copies")
seissol_add_comparison(
    NAME bitwise-box48-copies-states-x
    REFERENCE easi-box48-copies-x
    CANDIDATE easi-box48-copies-states
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_EVEN}
    LABELS "bitwise;copies")
seissol_add_comparison(
    NAME bitwise-box48-copies-states-y
    REFERENCE easi-box48-copies-y
    CANDIDATE easi-box48-copies-states
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_ODD}
    LABELS "bitwise;copies")

# An easi boundary condition per simulation: the Dirichlet faces of every
# second simulation take the other constant terms, in both copies of the
# cells and in the fused simulations of each, which share the map of their
# faces. The runs whose simulations all take the one or the other condition
# hold them alike.
set(SEISSOL_DIRICHLET ${SEISSOL_ZERO_STATE} "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}")
foreach (_copies_quantity v1 v2)
  seissol_add_minitest(
      NAME dirichlet-box48-copies-${_copies_quantity}
      MESH box48-dirichlet
      PARAMETERS wave.par
      STEPS 20
      SET ${SEISSOL_DIRICHLET} "equations.BoundaryFileName='boundary-${_copies_quantity}.yaml'"
      FILES boundary-${_copies_quantity}.yaml
      REQUIRE equations=elastic
      ACTIVITY "v1|v2:1e-12"
      NO_SNAPSHOT
      LABELS "smoke;copies")
endforeach()
set(SEISSOL_COPIES_BOUNDARIES ${SEISSOL_DIRICHLET}
    "equations.BoundaryFileName='boundary-v1.yaml'")
math(EXPR _copies_last "${SEISSOL_COPIES_SIMULATIONS} - 1")
foreach (_copies_simulation RANGE ${_copies_last})
  math(EXPR _copies_odd "${_copies_simulation} % 2")
  if (_copies_odd)
    list(APPEND SEISSOL_COPIES_BOUNDARIES
         "equations.BoundaryFileName${_copies_simulation}='boundary-v2.yaml'")
  endif()
endforeach()
seissol_add_minitest(
    NAME dirichlet-box48-copies-terms
    MESH box48-dirichlet
    PARAMETERS wave.par
    STEPS 20
    SET ${SEISSOL_COPIES_BOUNDARIES}
    FILES boundary-v1.yaml boundary-v2.yaml
    REQUIRE equations=elastic
    ACTIVITY "v1|v2:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;copies")
seissol_add_comparison(
    NAME bitwise-box48-copies-terms-v1
    REFERENCE dirichlet-box48-copies-v1
    CANDIDATE dirichlet-box48-copies-terms
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_EVEN}
    LABELS "bitwise;copies")
seissol_add_comparison(
    NAME bitwise-box48-copies-terms-v2
    REFERENCE dirichlet-box48-copies-v2
    CANDIDATE dirichlet-box48-copies-terms
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_ODD}
    LABELS "bitwise;copies")
# the faces of the simulations fused into a cell share the map of their
# condition, so their conditions cannot differ in it
seissol_add_minitest(
    NAME dirichlet-box48-fused-maps
    MESH box48-dirichlet
    PARAMETERS wave.par
    STEPS 20
    SET ${SEISSOL_ZERO_STATE} "equations.BoundaryFileName='boundary-v1.yaml'"
        "equations.BoundaryFileName1='boundary-mirror.yaml'"
    FILES boundary-v1.yaml boundary-mirror.yaml
    REQUIRE equations=elastic "fused>=2"
    EXPECT_FAILURE "can differ in their constant terms only"
    LABELS "smoke;copies")

# the simulations fused into a cell share its material, so they cannot name
# different files
seissol_add_minitest(
    NAME smoke-box48-fused-materials
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    SET "equations.MaterialFileName1='material-elastic-soft.yaml'"
    FILES material-elastic-soft.yaml
    REQUIRE equations=elastic "fused>=2"
    EXPECT_FAILURE "are fused into the same cell"
    LABELS "smoke;copies")

# A plasticity file per simulation: every second simulation takes parameters
# under which the wave yields, in both copies of the cells and in the fused
# simulations of each. The runs whose simulations all take the one or the other
# parameters hold them alike. The plastic strain relaxes within a step, so that
# yielding shows.
set(SEISSOL_COPIES_PLASTICITY "equations.Plasticity=1" "equations.Tv=6.103515625e-05"
    "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}")
seissol_add_minitest(
    NAME plasticity-skew48-copies
    MESH skew48
    PARAMETERS wave.par
    MATERIAL material-plastic.yaml
    STEPS 20
    SET ${SEISSOL_COPIES_PLASTICITY}
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    NO_SNAPSHOT
    LABELS "plasticity;copies")
# a simulation without a plasticity file of its own takes the one of the first simulation
set(SEISSOL_COPIES_WEAK ${SEISSOL_COPIES_PLASTICITY}
    "equations.PlasticityFileName='plasticity-weak.yaml'")
math(EXPR _copies_last "${SEISSOL_COPIES_SIMULATIONS} - 1")
seissol_add_minitest(
    NAME plasticity-skew48-copies-weak
    MESH skew48
    PARAMETERS wave.par
    MATERIAL material-plastic.yaml
    STEPS 20
    SET ${SEISSOL_COPIES_WEAK}
    FILES plasticity-weak.yaml
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    NO_SNAPSHOT
    LABELS "plasticity;copies")
set(SEISSOL_COPIES_PLASTICITIES ${SEISSOL_COPIES_PLASTICITY})
foreach (_copies_simulation RANGE ${_copies_last})
  math(EXPR _copies_odd "${_copies_simulation} % 2")
  if (_copies_odd)
    list(APPEND SEISSOL_COPIES_PLASTICITIES
         "equations.PlasticityFileName${_copies_simulation}='plasticity-weak.yaml'")
  endif()
endforeach()
seissol_add_minitest(
    NAME plasticity-skew48-copies-parameters
    MESH skew48
    PARAMETERS wave.par
    MATERIAL material-plastic.yaml
    STEPS 20
    SET ${SEISSOL_COPIES_PLASTICITIES}
    FILES plasticity-weak.yaml
    REQUIRE equations=elastic
    ACTIVITY "v1:1e-12"
    NO_SNAPSHOT
    LABELS "plasticity;copies")
seissol_add_comparison(
    NAME bitwise-skew48-copies-parameters-strong
    REFERENCE plasticity-skew48-copies
    CANDIDATE plasticity-skew48-copies-parameters
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_EVEN}
    LABELS "bitwise;plasticity;copies")
seissol_add_comparison(
    NAME bitwise-skew48-copies-parameters-weak
    REFERENCE plasticity-skew48-copies-weak
    CANDIDATE plasticity-skew48-copies-parameters
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_ODD}
    LABELS "bitwise;plasticity;copies")

# Point sources are dealt out to the simulations of a run in turn, and the copy
# of the cells that holds the simulation of a source takes it, with the density
# of its material. source-s{fused}.fsrm has a source for every fused
# simulation, source-copies-s{fused}.fsrm the same sources twice over (see
# sources.py): the simulations s and fused + s of the run of two copies carry
# the wave of the simulation s of the run of the fused simulations alone, which
# no other simulation carries. Where the second copy of the cells is denser,
# its simulations are those of the run whose copies all are.
set(SEISSOL_SOURCES ${SEISSOL_ZERO_STATE} "SourceType.Type=50")
set(SEISSOL_COPIES_SOURCES ${SEISSOL_SOURCES} "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}"
    "SourceType.FileName='source-copies-s{fused}.fsrm'")
seissol_add_minitest(
    NAME source-box48
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    SET ${SEISSOL_SOURCES} "SourceType.FileName='source-s{fused}.fsrm'"
    FILES source-s{fused}.fsrm
    REQUIRE equations=elastic
    ACTIVITY "v1|v2|v3:1e-12"
    LABELS "smoke;source")
seissol_add_minitest(
    NAME source-box48-copies
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    SET ${SEISSOL_COPIES_SOURCES}
    FILES source-copies-s{fused}.fsrm
    REQUIRE equations=elastic
    ACTIVITY "v1|v2|v3:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;source;copies")
seissol_add_minitest(
    NAME source-box48-copies-mpi4
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    RANKS 4
    SET ${SEISSOL_COPIES_SOURCES}
    FILES source-copies-s{fused}.fsrm
    REQUIRE equations=elastic
    ACTIVITY "v1|v2|v3:1e-12"
    LABELS "smoke;source;mpi;copies")
set(SEISSOL_SOURCES_SECOND "")
math(EXPR _copies_last "${SEISSOL_CAP_SIMULATIONS} - 1")
foreach (_copies_simulation RANGE ${_copies_last})
  math(EXPR _copies_candidate "${SEISSOL_CAP_SIMULATIONS} + ${_copies_simulation}")
  list(APPEND SEISSOL_SOURCES_SECOND "${_copies_simulation}:${_copies_candidate}")
endforeach()
seissol_add_comparison(
    NAME bitwise-box48-sources-first
    REFERENCE source-box48
    CANDIDATE source-box48-copies
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_FIRST}
    LABELS "bitwise;source;copies")
seissol_add_comparison(
    NAME bitwise-box48-sources-second
    REFERENCE source-box48
    CANDIDATE source-box48-copies
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_SOURCES_SECOND}
    LABELS "bitwise;source;copies")
seissol_add_comparison(
    NAME bitwise-box48-sources-mpi4
    REFERENCE source-box48-copies
    CANDIDATE source-box48-copies-mpi4
    TOLERANCE 0
    LABELS "bitwise;source;mpi;copies")
seissol_add_minitest(
    NAME source-box48-copies-dense
    MESH box48
    PARAMETERS wave.par
    MATERIAL material-elastic-dense.yaml
    STEPS 20
    SET ${SEISSOL_COPIES_SOURCES}
    FILES source-copies-s{fused}.fsrm
    REQUIRE equations=elastic
    ACTIVITY "v1|v2|v3:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;source;copies")
set(SEISSOL_SOURCES_MATERIALS ${SEISSOL_COPIES_SOURCES})
math(EXPR _copies_last "${SEISSOL_COPIES_SIMULATIONS} - 1")
foreach (_copies_simulation RANGE ${SEISSOL_CAP_SIMULATIONS} ${_copies_last})
  list(APPEND SEISSOL_SOURCES_MATERIALS
       "equations.MaterialFileName${_copies_simulation}='material-elastic-dense.yaml'")
endforeach()
seissol_add_minitest(
    NAME source-box48-copies-materials
    MESH box48
    PARAMETERS wave.par
    STEPS 20
    SET ${SEISSOL_SOURCES_MATERIALS}
    FILES source-copies-s{fused}.fsrm material-elastic-dense.yaml
    REQUIRE equations=elastic
    ACTIVITY "v1|v2|v3:1e-12"
    NO_SNAPSHOT
    LABELS "smoke;source;copies")
seissol_add_comparison(
    NAME bitwise-box48-sources-materials-first
    REFERENCE source-box48-copies
    CANDIDATE source-box48-copies-materials
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_FIRST}
    LABELS "bitwise;source;copies")
seissol_add_comparison(
    NAME bitwise-box48-sources-materials-second
    REFERENCE source-box48-copies-dense
    CANDIDATE source-box48-copies-materials
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_SECOND}
    LABELS "bitwise;source;copies")


# ---------------------------------------------------------------------------
# Several configurations in one run. ConfigMap puts the cells beyond the split
# plane of mat48 (group 2) into another configuration of the executable than
# the cells before it, so that each face of the plane joins the two. The pairs
# are those of the set ci-interaction (cmake/confignames.cmake): across their
# faces, the time integrals of a neighbor change order and precision, within a
# family (elastic, viscoelastic and anisotropic; acoustic and viscoacoustic) and
# between the families (a solid, a fluid and a poroelastic solid). A case runs in
# the configuration of the near half, in a build that holds both; the material
# file has a medium for each half. Only a
# pair in double precision keeps a snapshot reference: in another environment,
# a run that is partly single precision deviates from its reference by more
# than the snapshot tolerance of double precision allows.
# ---------------------------------------------------------------------------

foreach (pair
    elastic-o6-f64:acoustic-o4-f64:solid-fluid
    elastic-o6-f64:viscoelastic-m3-o5-f32:solid-solid
    elastic-o6-f64:anisotropic-o6-f64:solid-uniform
    anisotropic-o6-f64:viscoelastic-m3-o5-f32:solid-solid
    acoustic-o4-f64:elastic-o6-f64:fluid-solid
    acoustic-o4-f64:anisotropic-o6-f64:fluid-solid
    elastic-o6-f64:poroelastic-o5-f64:solid-poro
    poroelastic-o5-f64:acoustic-o4-f64:poro-fluid
    viscoelastic-m3-o5-f32:viscoacoustic-m3-o4-f32:solid-fluid
    acoustic-o4-f64:viscoacoustic-m3-o4-f32:fluid-fluid)
  string(REPLACE ":" ";" pair ${pair})
  list(GET pair 0 near)
  list(GET pair 1 far)
  list(GET pair 2 media)
  set(snapshot "")
  if (near MATCHES "-f32" OR far MATCHES "-f32")
    set(snapshot NO_SNAPSHOT)
  endif()
  set(interaction
      MESH mat48
      PARAMETERS wave.par
      MATERIAL material-interaction-${media}.yaml
      STEPS 20
      REQUIRE configuration=${near} configurations=${far}
      # the relaxation frequencies are read where a configuration of the run is anelastic
      SET "equations.ConfigMap='2:${far}'" "equations.FreqCentral=60.0" "equations.FreqRatio=100.0")
  seissol_add_minitest(
      NAME interaction-${near}-${far}
      ${interaction}
      ${snapshot}
      ACTIVITY "v1:1e-12"
      LABELS "smoke;interaction")
  seissol_add_minitest(
      NAME null-interaction-${near}-${far}
      ${interaction}
      SET ${SEISSOL_ZERO_STATE}
      EXPECT_ZERO
      LABELS "null;interaction")
endforeach()

# The solid of material-interaction-solid-uniform.yaml is isotropic, also as the stiffness tensor
# that an anisotropic configuration reads: with the far half anisotropic, the run solves the same
# problem as one with both halves elastic, apart from rounding, the Riemann problem at the plane
# included, which it poses in the anisotropic material. (A planar wave is set up from one cell per
# configuration, so that it is the same wave in both runs in a uniform medium only; the Riemann
# problem between different solids is held against the isotropic one by a unit test.) Compared
# sample by sample, the two runs differ by 1e-14 relative to the values, and by up to 2e-8 where a
# small stress component crosses zero; an error in the conversion or the Riemann problem gives a
# difference of order one.
seissol_add_minitest(
    NAME interaction-elastic-o6-f64-alone
    MESH mat48
    PARAMETERS wave.par
    MATERIAL material-interaction-solid-uniform.yaml
    STEPS 20
    REQUIRE configuration=elastic-o6-f64 configurations=anisotropic-o6-f64
    NO_SNAPSHOT
    ACTIVITY "v1:1e-12"
    LABELS "smoke;interaction")
seissol_add_comparison(
    NAME identity-interaction-anisotropic
    REFERENCE interaction-elastic-o6-f64-alone
    CANDIDATE interaction-elastic-o6-f64-anisotropic-o6-f64
    TOLERANCE 1e-6
    LABELS "identity;interaction")

# Likewise, the solid of material-interaction-fluid-solid.yaml is isotropic as a stiffness tensor:
# with the far half anisotropic, the run solves the same problem as with the far half elastic, the
# face to the fluid included, whose Riemann problem takes the impedance of the anisotropic solid
# from its Christoffel matrix.
seissol_add_comparison(
    NAME identity-interaction-anisotropic-fluid
    REFERENCE interaction-acoustic-o4-f64-elastic-o6-f64
    CANDIDATE interaction-acoustic-o4-f64-anisotropic-o6-f64
    TOLERANCE 1e-6
    LABELS "identity;interaction")

# Two poroelastic sandstones in the halves of mat48, in one configuration: with the pore pressure
# continuous across the boundary of their groups (open), with a skin at it (skin), and with a skin
# of a conductance so large (nearly-open) that it solves the same problem as open apart from
# rounding. The Riemann problem at a face with a skin is posed with the impedances of both sides,
# the one of open faces with the eigenvectors of the poroelastic material.
foreach (pores open nearly-open skin)
  set(poroelastic
      MESH mat48
      PARAMETERS wave.par
      MATERIAL material-interaction-poro-${pores}.yaml
      STEPS 20
      REQUIRE configuration=poroelastic-o5-f64)
  seissol_add_minitest(
      NAME interaction-poroelastic-o5-f64-${pores}
      ${poroelastic}
      NO_SNAPSHOT
      ACTIVITY "v1:1e-12"
      LABELS "smoke;interaction")
  seissol_add_minitest(
      NAME null-interaction-poroelastic-o5-f64-${pores}
      ${poroelastic}
      SET ${SEISSOL_ZERO_STATE}
      EXPECT_ZERO
      LABELS "null;interaction")
endforeach()
seissol_add_comparison(
    NAME identity-interaction-poroelastic-skin
    REFERENCE interaction-poroelastic-o5-f64-open
    CANDIDATE interaction-poroelastic-o5-f64-nearly-open
    TOLERANCE 1e-6
    LABELS "identity;interaction")

# ---------------------------------------------------------------------------
# Friction laws (seissol_add_rupture).
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

foreach (law aging:3 slip:4)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  seissol_add_rupture(NAME rupture-rs-${name} FL ${number} FAULT fault-rs.yaml
                      SET ${SEISSOL_RS})
  seissol_add_rupture(NAME rupture-rs-${name}-tp FL ${number} FAULT fault-rs-tp.yaml
                      SET ${SEISSOL_RS} ${SEISSOL_TP})
endforeach()
# a friction of 0.625 is beyond what this law holds at the initial slip rate in
# a non-negative state; see its fault files
seissol_add_rupture(NAME rupture-rs-svw FL 7 FAULT fault-rs-svw.yaml
                    SET ${SEISSOL_RS})
seissol_add_rupture(NAME rupture-rs-svw-tp FL 7 FAULT fault-rs-svw-tp.yaml
                    SET ${SEISSOL_RS} ${SEISSOL_TP})
seissol_add_rupture(NAME rupture-rs-fvw FL 103 FAULT fault-rs-fvw.yaml
                    SET ${SEISSOL_RS_FVW})
seissol_add_rupture(NAME rupture-rs-fvw-tp FL 103 FAULT fault-rs-fvw-tp.yaml
                    SET ${SEISSOL_RS_FVW} ${SEISSOL_TP})

foreach (law yoffe:33 gaussian:34)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  seissol_add_rupture(NAME rupture-imposed-${name} FL ${number} FAULT fault-imposed.yaml
                      SET "DynamicRupture.SlipRateOutputType=0")
endforeach()
# the delta pulse gets its onset off the step boundaries; see its fault file
seissol_add_rupture(NAME rupture-imposed-delta FL 35 FAULT fault-imposed-delta.yaml
                    SET "DynamicRupture.SlipRateOutputType=0")

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

# The energy output, which sums over the cells and the faces of a rank and then
# over the ranks, in another order on another number of them or of threads; the
# cases of its own keep the parameters of the others, and with them their
# snapshots. The momentum is reported, but not judged (see
# run_case.compare_energies).
set(SEISSOL_ENERGY_OUTPUT "Output.EnergyOutput=1" "Output.EnergyOutputInterval=9.765625e-04")
set(SEISSOL_ENERGY_SUMS 1e-13)
seissol_add_rupture(NAME rupture-lsw-energy FL 16 FAULT fault-lsw.yaml
                    SET ${SEISSOL_ENERGY_OUTPUT} NO_SNAPSHOT LABELS energy)
seissol_add_rupture(NAME rupture-lsw-energy-mpi4 FL 16 FAULT fault-lsw.yaml RANKS 4
                    SET ${SEISSOL_ENERGY_OUTPUT} LABELS "energy;mpi")
seissol_add_comparison(
    NAME identity-rupture-lsw-energy-mpi4
    REFERENCE rupture-lsw-energy
    CANDIDATE rupture-lsw-energy-mpi4
    TOLERANCE 0
    ENERGIES ${SEISSOL_ENERGY_SUMS}
    LABELS "identity;mpi;dr;energy")
seissol_add_rupture(NAME rupture-rs-fvw-tp-energy FL 103 FAULT fault-rs-fvw-tp.yaml
                    SET ${SEISSOL_RS_FVW} ${SEISSOL_TP} ${SEISSOL_ENERGY_OUTPUT}
                    NO_SNAPSHOT LABELS energy)
seissol_add_rupture(NAME rupture-rs-fvw-tp-energy-threads4 FL 103 FAULT fault-rs-fvw-tp.yaml
                    SET ${SEISSOL_RS_FVW} ${SEISSOL_TP} ${SEISSOL_ENERGY_OUTPUT} THREADS 4
                    LABELS "energy;threads")
seissol_add_comparison(
    NAME identity-rupture-rs-fvw-tp-energy-threads4
    REFERENCE rupture-rs-fvw-tp-energy
    CANDIDATE rupture-rs-fvw-tp-energy-threads4
    TOLERANCE 0
    ENERGIES ${SEISSOL_ENERGY_SUMS}
    LABELS "identity;threads;dr;energy")

# The clusters take their steps along the time stepping plan (SEISSOL_TIMESTEPPING_PLAN, see
# TimeManager::followPlan) instead of deciding them as actors: the same run, bit for bit, with the
# faces of a fault and through two epochs (the receivers write after 32 of the 64 steps).
seissol_add_rupture(NAME rupture-lsw-plan FL 16 FAULT fault-lsw.yaml NO_SNAPSHOT
                    SET "Output.ReceiverOutputInterval=0.001953125"
                    ENV "SEISSOL_TIMESTEPPING_PLAN=1"
                    LABELS plan)
seissol_add_comparison(
    NAME bitwise-rupture-lsw-plan
    REFERENCE rupture-lsw
    CANDIDATE rupture-lsw-plan
    TOLERANCE 0
    LABELS "bitwise;dr;plan")
# the faces of a layer are solved by the threads side by side, and a law that keeps the scratch of
# a face in itself shares it among them: thermal pressurization prepares coefficients for a face
# that the slip-rate solve reads afterwards
seissol_add_rupture(NAME rupture-rs-fvw-tp-threads4 FL 103 FAULT fault-rs-fvw-tp.yaml
                    SET ${SEISSOL_RS_FVW} ${SEISSOL_TP} THREADS 4 LABELS threads)
seissol_add_comparison(
    NAME bitwise-rupture-rs-fvw-tp-threads4
    REFERENCE rupture-rs-fvw-tp
    CANDIDATE rupture-rs-fvw-tp-threads4
    TOLERANCE 0
    LABELS "bitwise;threads;dr")

# a deprecated spelling of the same law has to give the same run
seissol_add_comparison(
    NAME identity-lsw-legacy
    REFERENCE rupture-lsw
    CANDIDATE rupture-lsw-legacy
    TOLERANCE 0
    LABELS "identity;dr")

# A nucleation takes its onset and its rise time per point from the fault file
# where it provides them, and from the parameter file otherwise. The same values
# have to give the same run either way. A face holds the onsets and rise times
# of all its stress sources together, which the fault file does not know of.
# The parameter file of the case reading them from the fault file has values
# that would show if they were taken.
seissol_add_rupture(NAME rupture-lsw-onset FL 16 FAULT fault-lsw.yaml
                    SET "DynamicRupture.s_0=2.44140625e-04" NO_SNAPSHOT)
seissol_add_rupture(NAME rupture-lsw-onset-easi FL 16 FAULT fault-lsw-onset.yaml
                    SET "DynamicRupture.s_0=1.0" "DynamicRupture.t_0=1.0" NO_SNAPSHOT)
seissol_add_comparison(
    NAME identity-lsw-onset-easi
    REFERENCE rupture-lsw-onset
    CANDIDATE rupture-lsw-onset-easi
    TOLERANCE 0
    LABELS "identity;dr")

# ---------------------------------------------------------------------------
# Stage 1: restart from a checkpoint. A restart in the middle of a rupture has
# to go on as the run does without one: the state of the fault and of the
# cells comes from the checkpoint. A restart continues the files of the volume
# receivers but starts those of the on-fault receivers anew, so the volume
# receivers, which see the rupture radiate, are compared.
# ---------------------------------------------------------------------------

seissol_add_rupture(NAME rupture-lsw-restart FL 16 FAULT fault-lsw.yaml RESTART 32
                    LABELS restart)
seissol_add_comparison(
    NAME bitwise-rupture-lsw-restart
    REFERENCE rupture-lsw
    CANDIDATE rupture-lsw-restart
    TOLERANCE 0
    RECEIVERS volume
    LABELS "bitwise;restart;dr")

# A checkpoint tells the copies of the cells and of the fault faces apart, so
# the restart of a run of copies goes on as that run does. It restarts only the
# simulations it was written with: in a run of other ones, the data of a copy
# would land on another.
seissol_add_rupture(NAME rupture-lsw-copies FL 16 FAULT fault-lsw.yaml
                    SET "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}" NO_SNAPSHOT
                    LABELS copies)
seissol_add_rupture(NAME rupture-lsw-copies-restart FL 16 FAULT fault-lsw.yaml
                    SET "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}" RESTART 32
                    LABELS "copies;restart")
seissol_add_comparison(
    NAME bitwise-rupture-lsw-copies-restart
    REFERENCE rupture-lsw-copies
    CANDIDATE rupture-lsw-copies-restart
    TOLERANCE 0
    RECEIVERS volume
    LABELS "bitwise;copies;restart;dr")
seissol_add_rupture(NAME rupture-lsw-copies-restart-fused FL 16 FAULT fault-lsw.yaml
                    SET "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}" RESTART 32
                    RESTART_SET "equations.SimCount=${SEISSOL_CAP_SIMULATIONS}"
                    EXPECT_FAILURE "restarts only the simulations it was written with"
                    LABELS "copies;restart")

# ---------------------------------------------------------------------------
# Stage 1: a fault per simulation. Every second simulation takes a fault that
# weakens over twice the slip, in both copies of the fault faces and in the
# fused simulations of each. The runs whose simulations all take the one or the
# other fault hold them alike, on the fault and off it, and so does the run on
# four ranks, where fault faces end up in the copy layer.
# ---------------------------------------------------------------------------

seissol_add_rupture(NAME rupture-lsw-copies-dc FL 16 FAULT fault-lsw-dc.yaml
                    SET "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}" NO_SNAPSHOT
                    LABELS copies)
set(SEISSOL_COPIES_FAULTS "equations.SimCount=${SEISSOL_COPIES_SIMULATIONS}")
math(EXPR _copies_last "${SEISSOL_COPIES_SIMULATIONS} - 1")
foreach (_copies_simulation RANGE ${_copies_last})
  math(EXPR _copies_odd "${_copies_simulation} % 2")
  if (_copies_odd)
    list(APPEND SEISSOL_COPIES_FAULTS
         "DynamicRupture.ModelFileName${_copies_simulation}='fault-lsw-dc.yaml'")
  endif()
endforeach()
seissol_add_rupture(NAME rupture-lsw-copies-faults FL 16 FAULT fault-lsw.yaml
                    FILES fault-lsw-dc.yaml SET ${SEISSOL_COPIES_FAULTS} NO_SNAPSHOT
                    LABELS copies)
seissol_add_rupture(NAME rupture-lsw-copies-faults-mpi4 FL 16 FAULT fault-lsw.yaml RANKS 4
                    FILES fault-lsw-dc.yaml SET ${SEISSOL_COPIES_FAULTS}
                    LABELS "copies;mpi")
seissol_add_comparison(
    NAME bitwise-rupture-lsw-copies-faults-lsw
    REFERENCE rupture-lsw-copies
    CANDIDATE rupture-lsw-copies-faults
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_EVEN}
    LABELS "bitwise;copies;dr")
seissol_add_comparison(
    NAME bitwise-rupture-lsw-copies-faults-dc
    REFERENCE rupture-lsw-copies-dc
    CANDIDATE rupture-lsw-copies-faults
    TOLERANCE 0
    SIMULATIONS ${SEISSOL_COPIES_ODD}
    LABELS "bitwise;copies;dr")
seissol_add_comparison(
    NAME bitwise-rupture-lsw-copies-faults-mpi4
    REFERENCE rupture-lsw-copies-faults
    CANDIDATE rupture-lsw-copies-faults-mpi4
    TOLERANCE 0
    LABELS "bitwise;mpi;copies;dr")

# ---------------------------------------------------------------------------
# Stage 2: friction laws against an exact solution (seissol_add_fault1d).
#
# On fault1d, a uniform fault stays uniform: nothing is transferred along it,
# each side radiates a plane shear wave, and every point of the fault solves the
# same ordinary differential equation, which fault1d.py integrates to far below
# the error of the case. The fault receivers are compared with it, and so are
# the receivers in the volume, which the wave reaches. What this cannot see is
# a rupture front: there is none on a fault that breaks everywhere at once.
#
# The inputs are those of the cases on fault48, slowed down where the friction
# law would change the slip rate within a few steps: a shear wave takes ten
# milliseconds through a cell of fault1d, and a slip rate that changes faster
# sends out a wave the mesh does not resolve. Like the convergence cases, these
# run for long, so they carry that label too.
#
# An anisotropic build runs them all with material-anisotropic.yaml, which is
# isotropic in effect, and adds media with tilted axes: there the radiation
# damping depends on the direction of the slip, which then leaves the direction
# of the traction, and the slip changes the normal stress.
# ---------------------------------------------------------------------------

set(SEISSOL_FAULT1D_TP_OUTPUT "Pickpoint.OutputMask=1 1 1 1 1 1 0 1 1 1 1 1")
# The elementwise output as well, every 64 steps and on faces split in four, which
# takes its values from the degrees of freedom at the end of a step, where the
# receivers take them from the expansion of the step in time; its deviations are
# those of the receivers, prefixed with elementwise-
set(SEISSOL_FAULT1D_ELEMENTWISE "DynamicRupture.OutputPointType=5"
    "Elementwise.printtimeinterval_sec=3.90625e-03" "Elementwise.refinement_strategy=2"
    "Elementwise.refinement=1")
set(SEISSOL_FAULT1D_ELEMENTWISE_MASK "Elementwise.OutputMask=1 1 1 1 1 1 0 1 1 1 1 0")
set(SEISSOL_FAULT1D_ELEMENTWISE_TP_MASK "Elementwise.OutputMask=1 1 1 1 1 1 0 1 1 1 1 1")

seissol_add_fault1d(NAME fault1d-lsw FL 16 FAULT fault1d-lsw.yaml
                    SET ${SEISSOL_FAULT1D_ELEMENTWISE} ${SEISSOL_FAULT1D_ELEMENTWISE_MASK}
                    LABELS convergence)
# the forced rupture weakens the fault over t_0
seissol_add_fault1d(NAME fault1d-lsw-forced FL 16 FAULT fault1d-lsw-forced.yaml
                    SET ${SEISSOL_NO_NUCLEATION} "DynamicRupture.t_0=3.125e-02"
                    LABELS convergence)
seissol_add_fault1d(NAME fault1d-lsw-tpapprox FL 1058 FAULT fault1d-lsw.yaml
                    LABELS convergence)
# the radiation damping of two media; across a fault that slips uniformly, the
# normal stress does not change, so the Prakash-Clifton regularization holds it
seissol_add_fault1d(NAME fault1d-lsw-bimaterial FL 6 FAULT fault1d-lsw.yaml
                    MESH fault1d-bimaterial
                    SET "DynamicRupture.pc_vstar=0.1" "DynamicRupture.pc_prakashLength=0.1"
                    LABELS convergence)
foreach (law aging:3 slip:4 svw:7)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  seissol_add_fault1d(NAME fault1d-rs-${name} FL ${number} FAULT fault1d-rs-${name}.yaml
                      STEPS 2048 SET ${SEISSOL_RS} LABELS convergence)
endforeach()
seissol_add_fault1d(NAME fault1d-rs-svw-tp FL 7 FAULT fault1d-rs-svw-tp.yaml
                    STEPS 2048 SET ${SEISSOL_RS} ${SEISSOL_TP} ${SEISSOL_FAULT1D_TP_OUTPUT}
                    LABELS convergence)
seissol_add_fault1d(NAME fault1d-rs-fvw FL 103 FAULT fault1d-rs-fvw.yaml
                    STEPS 2048 SET ${SEISSOL_RS_FVW} LABELS convergence)
seissol_add_fault1d(NAME fault1d-rs-fvw-tp FL 103 FAULT fault1d-rs-fvw-tp.yaml
                    STEPS 2048 SET ${SEISSOL_RS_FVW} ${SEISSOL_TP} ${SEISSOL_FAULT1D_TP_OUTPUT}
                        ${SEISSOL_FAULT1D_ELEMENTWISE} ${SEISSOL_FAULT1D_ELEMENTWISE_TP_MASK}
                    LABELS convergence)
# the initial pore pressure lowers the effective normal stress the initial state
# has to hold the traction under
seissol_add_fault1d(NAME fault1d-rs-fvw-tp-pressurized FL 103
                    FAULT fault1d-rs-fvw-tp-pressurized.yaml STEPS 2048
                    SET ${SEISSOL_RS_FVW} ${SEISSOL_TP} "DynamicRupture.tp_iniPressure=-80.0e6"
                        ${SEISSOL_FAULT1D_TP_OUTPUT}
                    LABELS convergence)
foreach (law yoffe:33 gaussian:34)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  seissol_add_fault1d(NAME fault1d-imposed-${name} FL ${number} FAULT fault1d-imposed.yaml
                      STEPS 2048 SET "DynamicRupture.SlipRateOutputType=0"
                      LABELS convergence)
endforeach()
seissol_add_fault1d(NAME fault1d-anisotropic-lsw FL 16 FAULT fault1d-lsw-oblique.yaml
                    MATERIAL material-fault1d-anisotropic.yaml REQUIRE equations=anisotropic
                    LABELS convergence)
seissol_add_fault1d(NAME fault1d-anisotropic-lsw-bimaterial FL 16
                    FAULT fault1d-lsw-oblique.yaml MESH fault1d-bimaterial
                    MATERIAL material-fault1d-anisotropic.yaml REQUIRE equations=anisotropic
                    LABELS convergence)
seissol_add_fault1d(NAME fault1d-anisotropic-rs-aging FL 3 FAULT fault1d-rs-aging-oblique.yaml
                    STEPS 2048 SET ${SEISSOL_RS}
                    MATERIAL material-fault1d-anisotropic.yaml REQUIRE equations=anisotropic
                    LABELS convergence)
# the receivers report the normal stress the imposed slip leaves; the shear
# tractions are those of the Godunov state, see docs/fault-output.rst
seissol_add_fault1d(NAME fault1d-anisotropic-imposed FL 33 FAULT fault1d-imposed-oblique.yaml
                    STEPS 2048 SET "DynamicRupture.SlipRateOutputType=0"
                    MATERIAL material-fault1d-anisotropic.yaml REQUIRE equations=anisotropic
                    LABELS convergence)

# the outputs of rupture1d.par and the rupture velocity Vr, which the fault
# receivers take from the gradient of the rupture times of the Gauss points of
# their face (see docs/fault-output.rst)
set(SEISSOL_RUPTURE_VELOCITY_OUTPUT "Pickpoint.OutputMask=1 1 1 1 1 1 1 1 1 1 1 0")

# A forced rupture front along y on fault2d, at 1.5 times the shear wave
# speed. A front that outruns the shear waves outruns all it sends out, and
# behind it, every point of the fault solves the 0D problem of fault1d,
# shifted to its rupture time, with the radiation damping of the antiplane
# plane waves inclined at the angle of the Mach cone (see fault1d.Front).
# This holds at a receiver until the waves of where the front starts or meets
# itself arrive, or the part of the inclined waves that the boundaries
# x = +-W/2 reflect; the comparison ends there. It takes the slip rate, the
# normal stress and the direction of the slip that the fault receivers
# compute at their position, from the jump of the velocity, and the volume
# receivers. The other quantities of the fault receivers come from the Gauss
# point of the fault nearest to them, a few metres away along the front,
# which reaches it up to 0.4 ms earlier or later; but for the rupture
# velocity of their face, which the front sets at its own speed wherever it
# is steady, so at all three receivers.
seissol_add_fault1d(NAME fault2d-supershear FL 16 FAULT fault2d-supershear.yaml
                    MESH fault2d STEPS 2048 FILES faultreceiver2d.dat receiver2d.dat
                    SET "DynamicRupture.t_0=2.0e-02" "DynamicRupture.SlipRateOutputType=0"
                        "Pickpoint.PPFileName='faultreceiver2d.dat'"
                        "Output.RFileName='receiver2d.dat'"
                        ${SEISSOL_RUPTURE_VELOCITY_OUTPUT}
                    LABELS convergence)
# with twelve times the cells of fault1d, the case runs for a few minutes;
# it gets four times the time the others get
math(EXPR SEISSOL_FAULT2D_TIMEOUT "4 * ${SEISSOL_MINITEST_TIMEOUT}")
set_tests_properties(verification.fault2d-supershear${SEISSOL_VERIFICATION_SUFFIX}
                     PROPERTIES TIMEOUT ${SEISSOL_FAULT2D_TIMEOUT})

# A spontaneous antiplane rupture with linear slip weakening on
# fault2d-periodic, which runs out from a nucleation to both sides along y,
# held against the spectral boundary integral method of antiplane.py. The mesh
# is periodic across the fault as well, which the method takes into account:
# from 77 ms on, what the fault sends out returns to it from the other side,
# and from about 137 ms on, the fronts meet their periodic images at
# y = +-400 m. The fault receivers lie on Gauss points of the fault, at the
# same place in their triangles every 100 m along y, from the nucleation to
# where the front arrives some 40 ms before the end, and the points of the
# method run through them. The comparison takes the slip rate, the normal
# stress and the direction of the slip that they compute at their position,
# what they take from their Gauss point, where they lie on one with the
# quadrature of the build (see run_case._on_gauss_points), and the velocity
# and the stress of the volume receivers, near the nucleation and along the
# path of the fronts, up to a third of the way to the next fault (see
# run_case.compare_antiplane).
set(SEISSOL_ANTIPLANE
    "DynamicRupture.t_0=2.0e-02" "DynamicRupture.SlipRateOutputType=0"
    "Pickpoint.PPFileName='faultreceiver2d-periodic.dat'"
    "Output.RFileName='receiver2d-periodic.dat'")
seissol_add_fault1d(NAME fault2d-antiplane FL 16 FAULT fault2d-antiplane.yaml
                    MESH fault2d-periodic STEPS 2560
                    FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                    SET ${SEISSOL_ANTIPLANE} ${SEISSOL_RUPTURE_VELOCITY_OUTPUT}
                    LABELS convergence)
# the comparison takes minutes as well, for the field at the volume receivers
set_tests_properties(verification.fault2d-antiplane${SEISSOL_VERIFICATION_SUFFIX}
                     verification.exact.fault2d-antiplane${SEISSOL_VERIFICATION_SUFFIX}
                     PROPERTIES TIMEOUT ${SEISSOL_FAULT2D_TIMEOUT})
# The same with each of the rate-and-state laws, under SEISSOL_RS and
# SEISSOL_RS_FVW, from their creep of 1e-16 m/s. Their fault files choose a
# and L, and for severe velocity weakening b, so that the friction rises by
# 9 to 18 MPa where the front arrives and weakens over a slip of centimetres,
# as in fault2d-antiplane.yaml, and the mesh resolves the fronts alike (see
# fault2d-antiplane-rs-aging.yaml).
foreach (law aging:3 slip:4 svw:7 fvw:103)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  if (number EQUAL 103)
    set(SEISSOL_ANTIPLANE_RS ${SEISSOL_RS_FVW})
  else()
    set(SEISSOL_ANTIPLANE_RS ${SEISSOL_RS})
  endif()
  seissol_add_fault1d(NAME fault2d-antiplane-rs-${name} FL ${number}
                      FAULT fault2d-antiplane-rs-${name}.yaml
                      MESH fault2d-periodic STEPS 2560
                      FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                      SET ${SEISSOL_ANTIPLANE_RS} ${SEISSOL_ANTIPLANE} LABELS convergence)
  set_tests_properties(verification.fault2d-antiplane-rs-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       verification.exact.fault2d-antiplane-rs-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       PROPERTIES TIMEOUT ${SEISSOL_FAULT2D_TIMEOUT})
endforeach()
# Fast velocity weakening with thermal pressurization, whose pore pressure
# rises by up to a fifth of the normal stress (see
# fault2d-antiplane-rs-fvw-tp.yaml); the fault receivers report it and the
# temperature of their Gauss points.
seissol_add_fault1d(NAME fault2d-antiplane-rs-fvw-tp FL 103
                    FAULT fault2d-antiplane-rs-fvw-tp.yaml MESH fault2d-periodic STEPS 2560
                    FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                    SET ${SEISSOL_RS_FVW} ${SEISSOL_TP} ${SEISSOL_FAULT1D_TP_OUTPUT}
                        ${SEISSOL_ANTIPLANE}
                    LABELS convergence)
set_tests_properties(verification.fault2d-antiplane-rs-fvw-tp${SEISSOL_VERIFICATION_SUFFIX}
                     verification.exact.fault2d-antiplane-rs-fvw-tp${SEISSOL_VERIFICATION_SUFFIX}
                     PROPERTIES TIMEOUT ${SEISSOL_FAULT2D_TIMEOUT})
# The Prakash-Clifton regularization of FL=6, which acts on its own here,
# between two halves of the same solid: the strength the slip rate sees
# follows that of linear slip weakening at the rate (|V| + v*) / L, which
# with L = 4 mm delays the front by 6 ms over its way to y = 350 m. And the
# proxy of thermal pressurization, with a fault file of its own.
seissol_add_fault1d(NAME fault2d-antiplane-lsw-bimaterial FL 6 FAULT fault2d-antiplane.yaml
                    MESH fault2d-periodic STEPS 2560
                    FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                    SET "DynamicRupture.pc_vstar=0.1" "DynamicRupture.pc_prakashLength=0.004"
                        ${SEISSOL_ANTIPLANE}
                    LABELS convergence)
seissol_add_fault1d(NAME fault2d-antiplane-lsw-tpapprox FL 1058
                    FAULT fault2d-antiplane-lsw-tpapprox.yaml MESH fault2d-periodic STEPS 2560
                    FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                    SET ${SEISSOL_ANTIPLANE} LABELS convergence)
# Instantaneous healing: once its slip rate has exceeded 1 m/s, the fault
# regains mu_s where it slows down below that, which locks it 5 to 8 ms
# before the stopping phases bring it to rest otherwise, with up to 5 % less
# slip. The run ends at 145 ms, 3 ms before the fault heals at the outer
# fault receivers, and 7 ms before the fronts that meet their images reload
# the healed fault there to its static strength. The fault receivers report
# the slip rate of the friction law (SlipRateOutputType=1), which the healed
# fault holds at zero, where the jump of the velocity across it rings after
# the sudden stop, by up to half a metre per second. The volume receivers are
# compared up to 91 ms, when the fault first heals (see
# run_case.compare_antiplane).
seissol_add_fault1d(NAME fault2d-antiplane-lsw-healing FL 16 FAULT fault2d-antiplane.yaml
                    MESH fault2d-periodic STEPS 2368
                    FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                    SET "DynamicRupture.lsw_healingThreshold=1.0" ${SEISSOL_ANTIPLANE}
                        "DynamicRupture.SlipRateOutputType=1"
                    LABELS convergence)
foreach (name bimaterial tpapprox healing)
  set_tests_properties(verification.fault2d-antiplane-lsw-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       verification.exact.fault2d-antiplane-lsw-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       PROPERTIES TIMEOUT ${SEISSOL_FAULT2D_TIMEOUT})
endforeach()
# An imposed slip, which sets in at y = 0 and runs out to both sides at
# 2.5 km/s (see fault2d-antiplane-imposed.yaml), at the rates of FL=33 and
# FL=34. FL=35 imposes all of the slip within the sub-step of the onset,
# which leaves a slip rate, and a field sent out, that depend on the time
# step and converge to no function, a delta in time; it has no comparison.
foreach (law yoffe:33 gaussian:34)
  string(REPLACE ":" ";" law ${law})
  list(GET law 0 name)
  list(GET law 1 number)
  seissol_add_fault1d(NAME fault2d-antiplane-imposed-${name} FL ${number}
                      FAULT fault2d-antiplane-imposed.yaml MESH fault2d-periodic STEPS 2560
                      FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                      SET ${SEISSOL_ANTIPLANE} LABELS convergence)
  set_tests_properties(verification.fault2d-antiplane-imposed-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       verification.exact.fault2d-antiplane-imposed-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       PROPERTIES TIMEOUT ${SEISSOL_FAULT2D_TIMEOUT})
endforeach()
# Spontaneous in-plane ruptures, of mode II, which slip along y, held against
# the in-plane spectral method of antiplane.py (antiplane.Inplane), whose P
# waves return across the period from 44 ms on: one that runs out below the
# shear wave speed and turns faster only where its fronts near their periodic
# images, and one that turns faster than the shear waves some 240 m from the
# nucleation (see fault2d-inplane.yaml and fault2d-inplane-supershear.yaml).
# The receivers are those of the antiplane cases.
foreach (name inplane inplane-supershear)
  seissol_add_fault1d(NAME fault2d-${name} FL 16 FAULT fault2d-${name}.yaml
                      MESH fault2d-periodic STEPS 2560
                      FILES faultreceiver2d-periodic.dat receiver2d-periodic.dat
                      SET ${SEISSOL_ANTIPLANE} LABELS convergence)
  set_tests_properties(verification.fault2d-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       verification.exact.fault2d-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       PROPERTIES TIMEOUT ${SEISSOL_FAULT2D_TIMEOUT})
endforeach()

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

# The planar wave in the tilted medium of the anisotropic fault1d cases, in
# which each of the three waves along x moves all components of the velocity.
# The analytical solution is an eigenvector of the Jacobian of the medium,
# whose rows of the stresses and of the velocities differ in scale by some
# thirteen orders of magnitude. Where the waves couple all components, it is
# accurate only if the eigensolver balances them; in the media of the cases
# above, they decouple. Its own thresholds, cases/thresholds-tilted.json,
# take the largest relative error over the quantities, which the small stress
# sigma_yz sets.
seissol_add_minitest(
    NAME convergence-periodic162-tilted
    MESH periodic162
    PARAMETERS wave.par
    MATERIAL material-fault1d-anisotropic.yaml
    STEPS 256
    THRESHOLDS thresholds-tilted.json
    REQUIRE equations=anisotropic
    NO_SNAPSHOT
    ACTIVITY "v1:1e-12"
    LABELS convergence)
