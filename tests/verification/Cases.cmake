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

seissol_add_fault1d(NAME fault1d-lsw FL 16 FAULT fault1d-lsw.yaml LABELS convergence)
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
# which reaches it up to 0.4 ms earlier or later.
seissol_add_fault1d(NAME fault2d-supershear FL 16 FAULT fault2d-supershear.yaml
                    MESH fault2d STEPS 2048 FILES faultreceiver2d.dat receiver2d.dat
                    SET "DynamicRupture.t_0=2.0e-02" "DynamicRupture.SlipRateOutputType=0"
                        "Pickpoint.PPFileName='faultreceiver2d.dat'"
                        "Output.RFileName='receiver2d.dat'"
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
                    SET ${SEISSOL_ANTIPLANE} LABELS convergence)
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
foreach (name bimaterial tpapprox)
  set_tests_properties(verification.fault2d-antiplane-lsw-${name}${SEISSOL_VERIFICATION_SUFFIX}
                       verification.exact.fault2d-antiplane-lsw-${name}${SEISSOL_VERIFICATION_SUFFIX}
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
