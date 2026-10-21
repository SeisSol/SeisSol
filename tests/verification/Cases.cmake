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
