# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Continuum damage-breakage rheology.

The state is a strain tensor, a velocity, and the two internal variables of the
rheology: the damage alpha and the breakage B. Neither internal variable takes
part in the face coupling -- they evolve through source terms local to the cell
-- so both sit in the layout without a face role.

The strain carries the traction role even though it is a strain rather than a
stress. The two transform alike under the face rotation, which is what the role
drives today; that the mechanical traction is a derived quantity here, obtained
from sigma(eps, alpha, B), is a distinction the face machinery does not yet
draw.

Stress is a convex blend of a solid and a granular branch, weighted by B. The
solid branch degrades with alpha through the coupling modulus gammaR; the
granular branch is a polynomial in the strain invariant ratio xi with four
coefficients. Damage and breakage grow with the strain energy once xi passes
the onset threshold, and the damage heals under compaction below it.
"""

import numpy as np
import yateto.functions as yf
from kernels.aderdg.nonlinearck import NonLinearCK
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind
from yateto import Scalar, Tensor
from yateto.memory import CSCMemoryLayout
from yateto.type import Datatype

#: Position of the two internal variables on the quantity axis.
ALPHA = 9
BREAKAGE = 10

#: Selects the trace of a symmetric tensor in Voigt order.
TRACE = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

#: Weights the Voigt components so a plain contraction gives eps_ij eps_ij.
VOIGT_SQUARE = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])


#: Rows of the directional flux fed by the velocity, as (source, target, factor).
#: The strain equations transport the symmetric velocity gradient, so the shear
#: rows pick up a half.
VELOCITY_FLUX = (
    ((0, 0, -1.0), (1, 3, -0.5), (2, 5, -0.5)),
    ((1, 1, -1.0), (0, 3, -0.5), (2, 4, -0.5)),
    ((2, 2, -1.0), (1, 4, -0.5), (0, 5, -0.5)),
)

#: Rows of the directional flux fed by the stress, as (Voigt source, target).
#: The momentum equations transport the traction of a face normal to the
#: direction, scaled by the inverse density.
STRESS_FLUX = (
    ((0, 6), (3, 7), (5, 8)),
    ((3, 6), (1, 7), (4, 8)),
    ((5, 6), (4, 7), (2, 8)),
)


def fluxMap(rows, sourceExtent, targetExtent, factor=-1.0):
    values = np.zeros((sourceExtent, targetExtent))
    for row in rows:
        source, target = row[0], row[1]
        values[source, target] = row[2] if len(row) > 2 else factor
    return values


#: The isotropic identity in Voigt, so that contracting it with a strain gives
#: that strain back: the shear rows carry a half because the Voigt pair counts
#: twice.
ISOTROPIC_VOIGT = np.diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])


def velocityRows(dim):
    """Velocity rows of one direction, in the quantity numbering."""
    return tuple((6 + row[0], row[1], row[2]) for row in VELOCITY_FLUX[dim])


def momentumMap():
    values = np.zeros((3, 6, 3))
    for dim in range(3):
        for voigt, velocity in STRESS_FLUX[dim]:
            values[dim, voigt, velocity - 6] = 1.0
    return values


def voigtLift(extent):
    values = np.zeros((6, extent))
    for i in range(6):
        values[i, i] = VOIGT_SQUARE[i]
    return values


def velocityLift(extent):
    values = np.zeros((3, extent))
    for i in range(3):
        values[i, 6 + i] = 1.0
    return values


def strainSelector(extent):
    values = np.zeros((extent, 6))
    for i in range(6):
        values[i, i] = 1.0
    return values


def unit(position, extent):
    values = np.zeros(extent)
    values[position] = 1.0
    return values


class DamageADERDG(NonLinearCK):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir)
        self.configure(matricesDir, memLayout, kwargs)

    #: Material parameters the kernels read per element, in the order they sit
    #: in the tensor that carries them.
    ParameterOrder = (
        "rhoInv",
        "lambda0",
        "mu0",
        "gammaR",
        "xi0",
        "damageRate",
        "breakageRate",
        "healingRate",
        "betaAlpha",
        "aB0",
        "aB1",
        "aB2",
        "aB3",
    )

    #: Both families are three matrices with one pattern each: the geometry a
    #: cell was meshed with, and the operator the recursion transports by.
    StarClones = {
        "star": ["star(0)", "star(1)", "star(2)"],
        "transport": ["transport(0)", "transport(1)", "transport(2)"],
    }

    def starMatrix(self, dim):
        """What the derivative recursion transports by.

        Not the geometry: the operator this material transports by follows the
        state, so it is assembled per step and the name that carries it is the
        assembled one.
        """
        return self.db.transport[dim]

    def primaryGroups(self):
        return [
            QuantityGroup("eps", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
            QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
            QuantityGroup("alpha", QuantityKind.INVARIANT),
            QuantityGroup("breakage", QuantityKind.INVARIANT),
        ]

    def name(self):
        return "damage"

    def addStepTensors(self):
        """Everything the step kernel reads or keeps between its statements.

        The nodal quantities are temporaries: one set, rewritten at every time
        node, so the kernel needs no scratch of its own beyond what yateto
        allocates for it. Only alpha and B survive a node, because they march.
        """
        nodes = self.num3DQuadraturePoints()
        nq = self.numQuantities()

        # The material is a property of the cell, and a scalar argument is
        # uniform over a batch, so the parameters arrive as one tensor per
        # element. Each is pulled out of it once per kernel with a constant
        # selector, into a value of rank zero; from there on it reads and
        # behaves exactly like a scalar.
        self.parameterNames = self.ParameterOrder
        count = len(self.parameterNames)
        self.materialParameters = Tensor("materialParameters", (count,))
        self.parameterPicks = {
            name: Tensor(
                f"pick{name[0].upper()}{name[1:]}", (count,), unit(index, count)
            )
            for index, name in enumerate(self.parameterNames)
        }
        self.parameterValues = {
            name: Tensor(name, (), temporary=True) for name in self.parameterNames
        }

        self.rhoInv = self.parameterValues["rhoInv"][""]
        self.lambda0 = self.parameterValues["lambda0"][""]
        self.mu0 = self.parameterValues["mu0"][""]
        self.gammaR = self.parameterValues["gammaR"][""]
        self.xi0 = self.parameterValues["xi0"][""]
        self.damageRate = self.parameterValues["damageRate"][""]
        self.breakageRate = self.parameterValues["breakageRate"][""]
        self.healingRate = self.parameterValues["healingRate"][""]
        self.betaAlpha = self.parameterValues["betaAlpha"][""]
        self.aB = [self.parameterValues[f"aB{i}"][""] for i in range(4)]

        # Guards xi against a vanishing second invariant. Tied to the working
        # precision rather than to the model, hence not a model parameter.
        self.floor = Scalar("invariantFloor")

        self.trace = Tensor("traceSelect", (6,), TRACE)
        self.voigt = Tensor("voigtSquare", (6,), VOIGT_SQUARE)
        self.pickAlpha = Tensor("selectAlpha", (nq,), unit(ALPHA, nq))
        self.pickBreakage = Tensor("selectBreakage", (nq,), unit(BREAKAGE, nq))
        self.epsInit = Tensor("epsInit", (6,))
        # The directional maps hold nine non-zeros at most, scattered, so a
        # bounding box buys nothing and they are stored by their pattern.
        self.toFluxV = [
            Tensor(
                f"velocityToFlux{axis}",
                (3, nq),
                fluxMap(VELOCITY_FLUX[d], 3, nq),
                CSCMemoryLayout,
            )
            for d, axis in enumerate("XYZ")
        ]
        self.toFluxS = [
            Tensor(
                f"stressToFlux{axis}",
                (6, nq),
                fluxMap(STRESS_FLUX[d], 6, nq),
                CSCMemoryLayout,
            )
            for d, axis in enumerate("XYZ")
        ]
        # Modale Koeffizienten auf den Zellmittelwert. Aus den Quadraturge-
        # wichten und der Auswertung an den Knoten gebildet, damit es nicht
        # davon abhaengt, wie die konstante Basisfunktion normiert ist.
        self.cellMean = Tensor(
            "cellMean", (self.num3DBasisFunctions(),), self.modalMeanWeights()
        )
        self.pickStrain = Tensor("selectStrain", (nq, 6), strainSelector(nq))
        # Wo die Geschwindigkeit eine Dehnungszeile speist, je Referenz-
        # richtung: reine Geometrie, ohne Modul.
        self.strainFlux = Tensor(
            "strainFlux",
            (3, nq, nq),
            np.stack([fluxMap(velocityRows(d), nq, nq) for d in range(3)]),
        )
        # Wo eine Spannungszeile eine Geschwindigkeit speist. Der Voigt-Index
        # der Zeile ist der des Paares (Geschwindigkeit, Richtung), was die
        # Tabelle der Spannungsfluesse bereits sagt.
        self.momentumPlace = Tensor("momentumPlace", (3, 6, 3), momentumMap())
        # Hebt den Geschwindigkeitsindex auf die Quantity-Achse. Erst hier,
        # damit die Zwischenergebnisse so schmal bleiben wie ihr Inhalt.
        self.liftVelocity = Tensor("liftVelocity", (3, nq), velocityLift(nq))
        self.deltaVoigt = Tensor("deltaVoigt", (6,), TRACE)
        self.isotropicVoigt = Tensor("isotropicVoigt", (6, 6), ISOTROPIC_VOIGT)
        # Traegt das Voigt-Gewicht und hebt im selben Zug den Dehnungsindex
        # auf die Quantity-Achse. Eine Konstante, damit das Muster mitkommt.
        self.voigtDiagonal = Tensor("voigtDiagonal", (6, nq), voigtLift(nq))
        self.unitColumn = Tensor("unitColumn", (1,), np.ones(1))
        # Picks the constant basis function, which is where a cell value goes
        # when it has to sit in a modal column.
        self.constantMode = Tensor(
            "constantMode",
            (self.num3DBasisFunctions(),),
            unit(0, self.num3DBasisFunctions()),
        )
        self.weights = Tensor("quadratureWeights", (nodes,), self.nodalMeanWeights())

        def temporary(name, columns=None, datatype=None):
            return self.nodalTensor(name, columns, temporary=True, datatype=datatype)

        self.eps = temporary("epsTotal", 6)
        self.i1 = temporary("invariantI1")
        self.i2 = temporary("invariantI2")
        self.rootI2 = temporary("rootI2")
        self.xi = temporary("xi")
        self.intact = temporary("intact")
        self.twoMuEff = temporary("twoMuEff")
        self.sigmaNodal = temporary("sigmaNodal", 6)
        self.critical = temporary("criticalDamage")
        self.drive = temporary("damageDrive")
        # A truth value, and declared as one: a temporary without a datatype
        # is the working precision, and then the condition written into it is
        # a bool where it is written and a double where it is read.
        self.growing = temporary("damageGrowing", datatype=Datatype.BOOL)
        # One per time node, not one reused: a node reaches its damage by
        # weighing every source sampled before it, so those have to still be
        # there. They are scalars per quadrature point, so the whole history
        # costs less than one nodal stress.
        self.sourceAlpha = [
            temporary(f"sourceAlpha{q}") for q in range(self.numTimeNodes())
        ]
        self.sourceBreakage = [
            temporary(f"sourceBreakage{q}") for q in range(self.numTimeNodes())
        ]
        # The integrals are accumulated at the nodes and projected once. The
        # projection is linear, so it commutes with the quadrature sum, and
        # doing it per node would run the same matrix over the same tensor
        # once per node for nothing.
        self.sigmaIntegral = temporary("sigmaIntegral", 6)
        # The stress of a node in modal form. Projected once and combined
        # afterwards: the projection is the expensive part, and every
        # coefficient of the expansion is a combination of the same ones.
        # One per variable, not two columns of one: a temporary written
        # through two subslices gets a buffer per statement, and then only the
        # last of them is the tensor.
        self.alphaModal = Tensor(
            "alphaModal", (self.num3DBasisFunctions(),), temporary=True
        )
        self.breakageModal = Tensor(
            "breakageModal", (self.num3DBasisFunctions(),), temporary=True
        )
        self.sigmaModal = OptionalDimTensor(
            "sigmaModal",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (self.num3DBasisFunctions(), 6),
            alignStride=True,
            temporary=True,
        )
        # One accumulator per source, not two columns of one. A temporary that
        # is written through two different subslices gets a buffer per
        # statement, and then only the last of them is the tensor: the two
        # columns would live in two places and the projection would read one.
        # The sources integrate into what the corrector lifts; the values
        # themselves integrate into what a face reads.
        self.alphaValueIntegral = temporary("alphaValueIntegral", 1)
        self.breakageValueIntegral = temporary("breakageValueIntegral", 1)
        self.alphaIntegral = temporary("alphaIntegral", 1)
        self.breakageIntegral = temporary("breakageIntegral", 1)
        self.meanAlpha = Tensor("meanAlpha", (1,), temporary=True)
        self.meanBreakage = Tensor("meanBreakage", (1,), temporary=True)

        # The two that march live across the nodes but not beyond the kernel.
        # The wave speed does: the face coupling is scaled with it.
        self.alphaNodal = temporary("alphaNodal")
        self.breakageNodal = temporary("breakageNodal")
        self.nodeWave = Tensor("nodeWaveSpeed", (1,), temporary=True)
        self.waveIntegral = Tensor("waveIntegral", (1,), temporary=True)
        self.nodeShear = Tensor("nodeShearSpeed", (1,), temporary=True)
        self.shearIntegral = Tensor("shearIntegral", (1,), temporary=True)

    def transportStrainOffset(self):
        return self.transportGroupSlice("eps")[0]

    def materialParameterNames(self):
        return self.parameterNames

    def parameterStatements(self):
        """Pulls every material parameter out of the cell's parameter tensor."""
        return [
            self.parameterValues[name][""]
            <= self.materialParameters["z"] * self.parameterPicks[name]["z"]
            for name in self.parameterNames
        ]

    def constitutiveStatements(self):
        """The material response at the nodes of :attr:`nodalState`.

        Everything that is a function of the state at one instant: the
        invariants, the effective moduli and the stress. Nothing marches and
        nothing is integrated, so the same block serves a timestep node and a
        single evaluation of the state.

        Expects `alphaNodal` and `breakageNodal` to hold the internal
        variables of the instant.
        """
        eps, i1, i2 = self.eps, self.i1, self.i2
        rootI2, xi, intact = self.rootI2, self.xi, self.intact
        alpha, breakage = self.alphaNodal, self.breakageNodal
        trace, floor = self.trace, self.floor
        mu0, lambda0, gammaR, xi0 = self.mu0, self.lambda0, self.gammaR, self.xi0
        aB = self.aB

        statements = []
        statements += [
            eps["lc"]
            <= yf.add(self.nodalState["lc"].subslice("c", 0, 6), self.epsInit["c"]),
            i1["l"] <= eps["lc"] * trace["c"],
            i2["l"] <= yf.mul(eps["lc"], eps["lc"]) * self.voigt["c"],
            rootI2["l"] <= yf.sqrt(yf.maximum(i2["l"], floor)),
            xi["l"]
            <= yf.where(
                yf.greater(i2["l"], floor),
                yf.div(i1["l"], rootI2["l"]),
                0.0,
            ),
            intact["l"] <= 1.0 - breakage["l"],
        ]

        # Stress: a convex blend of the solid and the granular branch. The
        # solid shear modulus degrades with alpha through gammaR; the granular
        # branch is the gradient of the potential
        #
        #     F = P(xi) I2,   P(xi) = aB0 + aB1 xi + aB2 xi^2 + aB3 xi^3,
        #
        # which is P'(xi) sqrt(I2) on the trace and 2 P(xi) - xi P'(xi) on the
        # strain. Written out with xi sqrt(I2) = I1, that is
        #
        #     (aB1 sqrt(I2) + 2 aB2 I1 + 3 aB3 xi I1) delta
        #       + (2 aB0 + aB1 xi - aB3 xi^3) eps,
        #
        # so the factors of aB0 and aB3 are the ones the differentiation puts
        # there and not the ones the polynomial carries.
        twoMuEff, sigma = self.twoMuEff, self.sigmaNodal
        statements += [
            twoMuEff["l"]
            <= 2.0 * mu0
            - 2.0 * gammaR * xi0 * alpha["l"]
            - gammaR * yf.mul(alpha["l"], xi["l"]),
            sigma["lc"] <= yf.mul(twoMuEff["l"], eps["lc"]),
            sigma["lc"]
            <= sigma["lc"]
            + (lambda0 * i1["l"] - gammaR * yf.mul(alpha["l"], rootI2["l"]))
            * trace["c"],
            sigma["lc"]
            <= yf.mul(intact["l"], sigma["lc"])
            + yf.mul(
                breakage["l"],
                yf.mul(
                    2.0 * aB[0]
                    + aB[1] * xi["l"]
                    - aB[3] * yf.mul(xi["l"], yf.mul(xi["l"], xi["l"])),
                    eps["lc"],
                ),
            ),
            sigma["lc"]
            <= sigma["lc"]
            + yf.mul(
                breakage["l"],
                2.0 * aB[2] * i1["l"]
                + 3.0 * aB[3] * yf.mul(xi["l"], i1["l"])
                + aB[1] * rootI2["l"],
            )
            * trace["c"],
        ]

        return statements

    def waveSpeedStatements(self):
        """The squares of the fastest wave of each family at this instant. The
        squares, because that is what the moduli are affine in; the root is
        taken once per face, by whoever reads a bound.

        Two of them, because a dissipation that scales every mode with the
        faster of the two spends the difference on the slower -- and a rupture
        is made of the slower ones. Both are speeds rather than impedances:
        each is formed with the density of the cell it belongs to, so a face
        may take the larger of what its two sides hand it without either
        needing the other's material.
        """
        # Floored, because a face takes a root of these. The effective shear
        # modulus is affine in the damage and falls through zero at the point
        # where the solid branch stops describing anything -- and the shear
        # square reaches it first, at twice the damage the compressional one
        # needs. Past it a root is a NaN, and a NaN in a bound is a NaN in
        # every flux that face carries.
        #
        # Zero is not a repair: a mode with no bound is a mode with no
        # dissipation, which is the wrong side to err on. It is what keeps a
        # cell that has left the model from taking the run with it, and the
        # state that got there is the thing to look at.
        pSquared = yf.maximum(self.rhoInv * (self.lambda0 + self.twoMuEff["l"]), 0.0)
        sSquared = yf.maximum(self.rhoInv * 0.5 * self.twoMuEff["l"], 0.0)
        return [
            self.nodeWave["u"] <= yf.mul(yf.max(pSquared, "l"), self.unitColumn["u"]),
            self.nodeShear["u"] <= yf.mul(yf.max(sSquared, "l"), self.unitColumn["u"]),
        ]

    def stepStatements(self, node, weight, width):
        """The material response at one time node, and what it contributes.

        ``weight`` is the quadrature weight of the node and ``width`` the
        width of the whole step; how far the internal variables have travelled
        by this node follows from the time nodes alone and arrives as
        literals.
        """
        first = node == 0
        i2 = self.i2
        xi, intact = self.xi, self.intact
        alpha, breakage = self.alphaNodal, self.breakageNodal
        floor = self.floor
        mu0, lambda0, gammaR, xi0 = self.mu0, self.lambda0, self.gammaR, self.xi0
        sigma = self.sigmaNodal

        statements = []

        if first:
            statements += self.parameterStatements()

        # The internal variables enter the step from the state and then march;
        # everything else is rebuilt at every node.
        if first:
            statements += [
                alpha["l"] <= self.nodalState["lp"] * self.pickAlpha["p"],
                breakage["l"] <= self.nodalState["lp"] * self.pickBreakage["p"],
            ]
        else:
            # The distance from the node before this one, which is the
            # difference of two rows of the rule. Written as a step from the
            # node before rather than as a reach from the start of the step:
            # an addition onto the variable that travels is one accumulation
            # into the place it lives, where a sum that begins somewhere else
            # would accumulate into that somewhere else instead.
            travel = self.timeMarch()
            step = travel[node] - travel[node - 1]
            for q in range(node):
                if step[q] == 0.0:
                    continue
                scale = float(step[q]) * width
                statements += [
                    alpha["l"] <= alpha["l"] + scale * self.sourceAlpha[q]["l"],
                    breakage["l"]
                    <= breakage["l"] + scale * self.sourceBreakage[q]["l"],
                ]

        statements += self.constitutiveStatements()

        # The cell as a whole: the means the source guard asks for, and the
        # square of its fastest wave at this node. The square, because that is
        # what the moduli are affine in; the root is taken once per face, by
        # whoever reads the integral.
        statements += [
            self.meanAlpha["u"]
            <= alpha["l"] * self.weights["l"] * self.unitColumn["u"],
            self.meanBreakage["u"]
            <= breakage["l"] * self.weights["l"] * self.unitColumn["u"],
        ]
        statements += self.waveSpeedStatements()
        statements += [
            (
                accumulator["u"] <= node["u"]
                if first
                else accumulator["u"] <= yf.maximum(accumulator["u"], node["u"])
            )
            for accumulator, node in (
                (self.waveIntegral, self.nodeWave),
                (self.shearIntegral, self.nodeShear),
            )
        ]

        # The critical damage at which breakage sets in: the smaller root of a
        # quadratic in alpha, capped against the modulus ratio and against one.
        quadA = (
            3.0 * gammaR * gammaR * yf.mul(xi["l"], xi["l"])
            - 3.0 * gammaR * gammaR
            + 6.0 * gammaR * gammaR * xi0 * xi["l"]
            + 4.0 * gammaR * gammaR * xi0 * xi0
        )
        quadB = (
            -(8.0 * mu0 + 6.0 * lambda0) * gammaR * xi0
            - gammaR * lambda0 * yf.mul(xi["l"], yf.mul(xi["l"], xi["l"]))
            - 6.0 * gammaR * mu0 * xi["l"]
        )
        quadC = 4.0 * mu0 * mu0 + 6.0 * mu0 * lambda0
        fromRoot = yf.div(
            -quadB
            - yf.sqrt(yf.maximum(yf.mul(quadB, quadB) - 4.0 * quadA * quadC, floor)),
            2.0 * quadA,
        )
        fromModuli = yf.div(2.0 * mu0, gammaR * (xi["l"] + 2.0 * xi0))

        critical = self.critical
        statements += [
            critical["l"]
            <= yf.minimum(
                yf.minimum(
                    yf.where(yf.greater(fromRoot, floor), fromRoot, 1.0),
                    yf.where(yf.greater(fromModuli, floor), fromModuli, 1.0),
                ),
                1.0,
            ),
        ]

        # Damage and breakage grow with the strain energy above the onset
        # threshold; below it the damage heals at its own rate, and the cell
        # stops once either variable saturates.
        drive, growing = self.drive, self.growing
        statements += [
            drive["l"] <= gammaR * yf.mul(intact["l"], yf.mul(i2["l"], xi["l"] + xi0)),
            growing["l"]
            <= yf.logical_and(
                yf.greater(xi["l"] + xi0, floor),
                yf.logical_and(
                    yf.less(yf.sum(self.meanAlpha["u"], "u"), 1.0),
                    yf.less(yf.sum(self.meanBreakage["u"], "u"), 1.0),
                ),
            ),
            self.sourceAlpha[node]["l"]
            <= yf.where(
                growing["l"],
                self.damageRate * drive["l"],
                self.healingRate * drive["l"],
            ),
            # Selected rather than multiplied by. `growing` is a truth value,
            # and a truth value that is also a factor is a tensor with two
            # datatypes: the host path promotes it and the device path derives
            # one type per occurrence, so the two occurrences disagree and the
            # generator refuses the kernel. Every reader of it is a condition
            # now, as the one above already was.
            self.sourceBreakage[node]["l"]
            <= yf.where(
                growing["l"],
                self.breakageRate
                * yf.mul(
                    yf.div(
                        1.0,
                        1.0
                        + yf.exp(yf.div(critical["l"] - alpha["l"], self.betaAlpha)),
                    ),
                    drive["l"],
                ),
                0.0,
            ),
        ]

        # What leaves the node: the stress and the source under the quadrature
        # weight, and the two internal variables marched to the next node.
        accSigma = self.sigmaIntegral
        projection = self.timeProjection()
        statements += [
            (
                accSigma["lc"] <= weight * sigma["lc"]
                if first
                else accSigma["lc"] <= accSigma["lc"] + weight * sigma["lc"]
            )
        ]
        for accumulator, source in (
            (self.alphaIntegral, self.sourceAlpha[node]),
            (self.breakageIntegral, self.sourceBreakage[node]),
            (self.alphaValueIntegral, alpha),
            (self.breakageValueIntegral, breakage),
        ):
            statements += [
                (
                    accumulator["ln"] <= weight * source["l"] * self.unitColumn["n"]
                    if first
                    else accumulator["ln"]
                    <= accumulator["ln"] + weight * source["l"] * self.unitColumn["n"]
                )
            ]
        # The expansion of what the cell carries beyond its state, one
        # weighted sum per coefficient. The weights are the projection onto
        # the Legendre basis and constant, so they arrive as literals.
        modal = self.db.projectQP[self.t("kl")]
        statements += [
            self.sigmaModal["kc"] <= modal * sigma["lc"],
            # The two internal variables ride along, because the moduli a wave
            # sees at a face follow them and a face cannot see the state.
            self.alphaModal["k"] <= modal * alpha["l"],
            self.breakageModal["k"] <= modal * breakage["l"],
        ]

        # Relative to the carried block, which begins where the columns shared
        # with the state end.
        base = self.transportStateExtent()

        def carried(name):
            start, stop = self.transportGroupSlice(name)
            return (start - base, stop - base)

        stress = carried("sigma")
        internalColumns = (
            (carried("alpha"), self.alphaModal),
            (carried("breakage"), self.breakageModal),
        )
        for i in range(self.order):
            # a numpy scalar is not a Python float, and yateto takes the latter
            weightOf = float(projection[i, node])
            if weightOf == 0.0:
                # A node this coefficient does not weigh contributes nothing to
                # it, and the way to say so is to leave the statement out. A
                # term scaled by zero is one the generator is entitled to drop,
                # and a dropped statement that was going to write leaves its
                # target partly unwritten -- which the generator then refuses,
                # having been told to write a region it no longer covers.
                continue
            # Whichever node this coefficient is first weighed at opens it; the
            # ones after add to it. Not the first node as such, because that
            # one may be a node the coefficient does not weigh.
            opens = not any(projection[i, q] != 0.0 for q in range(node))
            target = self.transportDer[i]
            statements += [
                (
                    target["kc"].subslice("c", *stress)
                    <= weightOf * self.sigmaModal["kc"]
                    if opens
                    else target["kc"].subslice("c", *stress)
                    <= target["kc"].subslice("c", *stress)
                    + weightOf * self.sigmaModal["kc"]
                )
            ]
            for columns, value in internalColumns:
                statements += [
                    (
                        target["kc"].subslice("c", *columns)
                        <= weightOf * value["k"] * self.unitColumn["c"]
                        if opens
                        else target["kc"].subslice("c", *columns)
                        <= target["kc"].subslice("c", *columns)
                        + weightOf * value["k"] * self.unitColumn["c"]
                    )
                ]

        return statements

    def finishStatements(self):
        """What the step leaves behind, projected back to modal form once."""
        stress = self.transportGroupSlice("sigma")
        shared = self.transportStateExtent()
        bounds = (
            (self.transportGroupSlice("waveIntegral"), self.waveIntegral),
            (self.transportGroupSlice("shearIntegral"), self.shearIntegral),
        )
        projection = self.db.projectQP[self.t("kl")]
        return (
            [
                self.I["kc"].subslice("c", *stress)
                <= projection * self.sigmaIntegral["lc"],
                self.I["kc"].subslice("c", *self.transportGroupSlice("alpha"))
                <= projection * self.alphaValueIntegral["lc"],
                self.I["kc"].subslice("c", *self.transportGroupSlice("breakage"))
                <= projection * self.breakageValueIntegral["lc"],
                # The two columns meet on the way out, where the target is a real
                # buffer and a subslice of it is a place rather than a binding.
                self.sourceI["kn"].subslice("n", 0, 1)
                <= projection * self.alphaIntegral["ln"],
                self.sourceI["kn"].subslice("n", 1, 2)
                <= projection * self.breakageIntegral["ln"],
            ]
            + [
                # A cell value in a modal column is the coefficient of the
                # constant basis function and nothing else.
                self.I["kc"].subslice("c", *columns)
                <= self.constantMode["k"] * value["c"]
                for columns, value in bounds
            ]
            + [
                # The bound is one number for the whole step, not a function of
                # time within it: a maximum does not restrict to a subinterval the
                # way an integral does. So it stands in the constant coefficient
                # of the carried expansion. The rest have no column for one --
                # their layout says so, rather than a statement writing a zero
                # into it -- so a neighbour reconstructing a subinterval reads
                # the step's bound, which overestimates, which is the side
                # Rusanov may err on.
                self.transportDer[0]["kc"].subslice(
                    "c", columns[0] - shared, columns[1] - shared
                )
                <= self.constantMode["k"] * value["c"]
                for columns, value in bounds
            ]
        )

    def addTransport(self, generator, targets):
        """The operator the derivative recursion transports by, per step.

        The moduli of this material follow the state, so what the recursion
        carries is a linearisation about where the cell currently is. It is
        taken at the cell mean, which is one modal coefficient away, and it is
        the tangent of the stress rather than its ratio to the strain -- a
        wave travels at the former.

        With n the mean strain normalised in the Frobenius norm, xi = tr(n)
        and g = gammaR alpha, the tangent is

          C = lambda0 d(x)d + (2 mu0 - 2 g xi0 - g xi) Isym
                - g (d(x)n + n(x)d) + g xi n(x)n,

        which is not isotropic: the last two groups have no pair of Lame
        parameters behind them. Written in Voigt it is a plain 6 by 6, and the
        directional operator is that matrix placed in the rows where a stress
        feeds a velocity, plus the geometry where a velocity feeds a strain.
        The Jacobian of the cell weighs the three reference directions.
        """
        nq = self.numQuantities()
        mean = Tensor("meanState", (nq,), temporary=True)
        meanStrain = Tensor("meanStrain", (6,), temporary=True)
        meanAlpha = Tensor("meanAlpha", (), temporary=True)
        direction = Tensor("strainDirection", (6,), temporary=True)
        tangent = Tensor("tangent", (6, 6), temporary=True)
        invariantI1 = Tensor("meanI1", (), temporary=True)
        invariantI2 = Tensor("meanI2", (), temporary=True)
        ratio = Tensor("meanXi", (), temporary=True)
        coupling = Tensor("meanCoupling", (), temporary=True)
        shear = Tensor("meanShear", (), temporary=True)
        place = [Tensor(f"momentumRows({r})", (6, 3), temporary=True) for r in range(3)]

        delta, isotropic = self.deltaVoigt, self.isotropicVoigt
        floor = self.floor

        statements = self.parameterStatements()
        statements += [
            mean["p"] <= self.cellMean["k"] * self.Q["kp"],
            meanStrain["c"] <= mean["p"] * self.pickStrain["pc"] + self.epsInit["c"],
            meanAlpha[""] <= mean["p"] * self.pickAlpha["p"],
            invariantI1[""] <= meanStrain["c"] * self.trace["c"],
            invariantI2[""]
            <= yf.mul(meanStrain["c"], meanStrain["c"]) * self.voigt["c"],
            # The direction of the strain, which is all the tangent asks of it.
            # A cell at rest has no direction to give, and there the material
            # is the undamaged one whatever the damage says.
            ratio[""]
            <= yf.where(
                yf.greater(invariantI2[""], floor),
                invariantI1[""] / yf.sqrt(invariantI2[""]),
                0.0,
            ),
            direction["c"]
            <= yf.where(
                yf.greater(invariantI2[""], floor),
                meanStrain["c"] / yf.sqrt(yf.maximum(invariantI2[""], floor)),
                0.0,
            ),
            coupling[""]
            <= yf.where(
                yf.greater(invariantI2[""], floor),
                self.gammaR * meanAlpha[""],
                0.0,
            ),
            shear[""]
            <= 2.0 * self.mu0
            - 2.0 * yf.mul(coupling[""], self.xi0)
            - yf.mul(coupling[""], ratio[""]),
            tangent["cd"]
            <= self.lambda0 * delta["c"] * delta["d"]
            + yf.mul(shear[""], isotropic["cd"])
            - yf.mul(coupling[""], delta["c"] * direction["d"])
            - yf.mul(coupling[""], direction["c"] * delta["d"])
            + yf.mul(yf.mul(coupling[""], ratio[""]), direction["c"] * direction["d"]),
        ]
        for r in range(3):
            statements += [
                place[r]["ci"] <= self.db.star[r]["d"] * self.momentumPlace["dci"],
                self.db.transport[r]["qp"]
                <= self.db.star[r]["d"] * self.strainFlux["dqp"]
                # The Voigt pair of a shear row counts twice in the
                # contraction the recursion performs, so the weight rides
                # along; both lifts are constants, which is what keeps this
                # inside the pattern the operator is stored by.
                - self.rhoInv
                * tangent["ce"]
                * self.voigtDiagonal["eq"]
                * place[r]["ci"]
                * self.liftVelocity["ip"],
            ]

        for target in targets:
            prefix = generate_kernel_name_prefix(target)
            generator.add(f"{prefix}damageTransport", statements, target=target)

    def addStateToTransport(self, generator, targets):
        """What a damaged cell transports, at one instant, from its state.

        The strain, the velocity and the two internal variables are the
        state's own columns, moved into the places the transport layout keeps
        them. The stress is not: it is a function of the state, so it is
        evaluated at the nodes and projected back -- an interpolation rather
        than a projection, because a nonlinear function of a polynomial is
        not one, which is the approximation the step makes as well.

        The bound a face scales its dissipation with is the speed of this
        instant, which is what a step would leave behind if it were the only
        node in it.
        """
        coupled = (0, self.transportStateExtent())
        stress = self.transportGroupSlice("sigma")
        wave = self.transportGroupSlice("waveIntegral")
        projection = self.db.projectQP[self.t("kl")]

        statements = self.parameterStatements()
        statements += [
            self.nodalState["lp"] <= self.db.evalAtQP[self.t("lk")] * self.Q["kp"],
            self.alphaNodal["l"] <= self.nodalState["lp"] * self.pickAlpha["p"],
            self.breakageNodal["l"] <= self.nodalState["lp"] * self.pickBreakage["p"],
        ]
        statements += self.constitutiveStatements()
        statements += self.waveSpeedStatements()
        statements += [
            self.sigmaModal["kc"] <= projection * self.sigmaNodal["lc"],
            self.I["kc"].subslice("c", *coupled)
            <= self.Q["kc"].subslice("c", *coupled),
            self.I["kc"].subslice("c", *stress) <= self.sigmaModal["kc"],
            self.I["kc"].subslice("c", *self.transportGroupSlice("alpha"))
            <= self.Q["kc"].subslice("c", ALPHA, ALPHA + 1),
            self.I["kc"].subslice("c", *self.transportGroupSlice("breakage"))
            <= self.Q["kc"].subslice("c", BREAKAGE, BREAKAGE + 1),
            self.I["kc"].subslice("c", *wave)
            <= self.constantMode["k"] * self.nodeWave["c"],
            self.I["kc"].subslice("c", *self.transportGroupSlice("shearIntegral"))
            <= self.constantMode["k"] * self.nodeShear["c"],
        ]

        for target in targets:
            prefix = generate_kernel_name_prefix(target)
            generator.add(f"{prefix}stateToTransport", statements, target=target)

    def addCellIntegral(self, generator, target, prefix):
        """What the cell adds to its own state, from its own integrals.

        The flux is linear in the velocity and the stress, and both are
        transported, so the volume term is a constant map: no nodal detour and
        nothing to evaluate. The source integral of the internal variables
        rides along, because it is the other half of the same update and it is
        addition.
        """
        velocity = self.transportGroupSlice("v")
        stress = self.transportGroupSlice("sigma")
        internal = (self.transportStateExtent(), self.numQuantities())

        volume = self.Q["kp"]
        for d in range(3):
            volume = volume + self.db.kDivM[d][self.t("kl")] * (
                self.I["lm"].subslice("m", *velocity) * self.toFluxV[d]["mp"]
                + self.rhoInv
                * self.I["lc"].subslice("c", *stress)
                * self.toFluxS[d]["cp"]
            )

        generator.add(
            f"{prefix}damageCellIntegral",
            [
                *self.parameterStatements(),
                self.Q["kp"] <= volume,
                self.Q["kn"].subslice("n", *internal)
                <= self.Q["kn"].subslice("n", *internal) + self.sourceI["kn"],
            ],
            target=target,
        )

    def fluxPattern(self):
        """Which transported quantity feeds which equation across a face.

        The strain equations transport the velocity, the momentum equations
        the stress, and the two internal variables neither. Those are the same
        two tables the volume term uses, read with the transport layout's
        offsets.
        """
        velocity = self.transportGroupSlice("v")[0]
        stress = self.transportGroupSlice("sigma")[0]
        pairs = []
        for direction in range(3):
            for row in VELOCITY_FLUX[direction]:
                pairs.append((velocity + row[0], row[1]))
            for row in STRESS_FLUX[direction]:
                pairs.append((stress + row[0], row[1]))
        return pairs

    def drFluxSolverStatements(self, fluxScale, fluxSolver):
        """A fault's flux solver: the flux of what the friction imposed.

        The imposed state arrives in the frame of the face, so the normal is
        the first axis and the tables reduce to their normal slice; what comes
        out is rotated back by the quantity rotation, as the star-matrix
        version does. The rows it lands on are the ones the flux of this system
        has: an imposed traction enters the momentum of the cell, an imposed
        velocity enters its strain -- which is why this cannot be the flux of a
        state, where the traction is a state variable and lands on itself.
        """
        nq = self.numQuantities()
        velocityMap = Tensor(
            "velocityFluxMap",
            (3, 3, nq),
            np.stack([fluxMap(VELOCITY_FLUX[d], 3, nq) for d in range(3)]),
        )
        stressMap = Tensor(
            "stressFluxMap",
            (3, 6, nq),
            np.stack([fluxMap(STRESS_FLUX[d], 6, nq) for d in range(3)]),
        )
        faceNormal = Tensor("faceFrameNormal", (3,), np.array([1.0, 0.0, 0.0]))
        rhoInv = Scalar("rhoInv")
        velocity = self.transportGroupSlice("v")
        stress = self.transportGroupSlice("sigma")

        return [
            fluxSolver["qp"].subslice("q", *velocity)
            <= fluxScale * velocityMap["dqk"] * faceNormal["d"] * self.T["pk"],
            fluxSolver["qp"].subslice("q", *stress)
            <= fluxScale * rhoInv * stressMap["dqk"] * faceNormal["d"] * self.T["pk"],
        ]

    def fluxSolverStatements(self, fluxScale, normal):
        """The constant half of the flux solver: the flux of the face normal.

        Half the flux of each side, which is the average the numerical flux
        takes; the jump that goes with it is the dissipation, and it is added
        per timestep. The velocity block is the material's only by the inverse
        density, which the stress block carries.
        """
        nq = self.numQuantities()
        velocityMap = Tensor(
            "velocityFluxMap",
            (3, 3, nq),
            np.stack([fluxMap(VELOCITY_FLUX[d], 3, nq) for d in range(3)]),
        )
        stressMap = Tensor(
            "stressFluxMap",
            (3, 6, nq),
            np.stack([fluxMap(STRESS_FLUX[d], 6, nq) for d in range(3)]),
        )
        rhoInv = Scalar("rhoInv")
        velocity = self.transportGroupSlice("v")
        stress = self.transportGroupSlice("sigma")

        return [
            self.fluxConstant["qp"].subslice("q", *velocity)
            <= 0.5 * fluxScale * velocityMap["dqp"] * normal["d"],
            self.fluxConstant["qp"].subslice("q", *stress)
            <= 0.5 * fluxScale * rhoInv * stressMap["dqp"] * normal["d"],
        ]


def kernel_class(**kwargs):
    solver = kwargs["solver"].lower()
    if solver == "nonlinearck":
        return DamageADERDG
    raise NotImplementedError(f"{solver} cannot advance a material with damage.")
