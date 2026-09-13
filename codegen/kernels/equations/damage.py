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
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind
from yateto import Scalar, Tensor
from yateto.memory import CSCMemoryLayout

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
        self.unitColumn = Tensor("unitColumn", (1,), np.ones(1))
        # Picks the constant basis function, which is where a cell value goes
        # when it has to sit in a modal column.
        self.constantMode = Tensor(
            "constantMode",
            (self.num3DBasisFunctions(),),
            unit(0, self.num3DBasisFunctions()),
        )
        self.weights = Tensor("quadratureWeights", (nodes,), self.nodalMeanWeights())

        def temporary(name, columns=None):
            return self.nodalTensor(name, columns, temporary=True)

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
        self.growing = temporary("damageGrowing")
        self.sourceAlpha = temporary("sourceAlpha")
        self.sourceBreakage = temporary("sourceBreakage")
        # The integrals are accumulated at the nodes and projected once. The
        # projection is linear, so it commutes with the quadrature sum, and
        # doing it per node would run the same matrix over the same tensor
        # once per node for nothing.
        self.sigmaIntegral = temporary("sigmaIntegral", 6)
        self.sourceIntegral = temporary("sourceIntegral", 2)
        self.meanAlpha = Tensor("meanAlpha", (1,), temporary=True)
        self.meanBreakage = Tensor("meanBreakage", (1,), temporary=True)

        # The two that march live across the nodes but not beyond the kernel.
        # The wave speed does: the face coupling is scaled with it.
        self.alphaNodal = temporary("alphaNodal")
        self.breakageNodal = temporary("breakageNodal")
        self.waveIntegral = Tensor("waveIntegral", (1,), temporary=True)
        self.intervalLength = Tensor("intervalLength", (1,), temporary=True)

    def materialParameterNames(self):
        return self.parameterNames

    def parameterStatements(self):
        """Pulls every material parameter out of the cell's parameter tensor."""
        return [
            self.parameterValues[name][""]
            <= self.materialParameters["z"] * self.parameterPicks[name]["z"]
            for name in self.parameterNames
        ]

    def stepStatements(self, node, weight, march):
        """The material response at one time node, and what it contributes.

        ``weight`` is the quadrature weight of the node, ``march`` the step
        from it to the next one. Both are scalars the caller sets, so which
        rule is being used is a property of the launch code and not of the
        kernel.
        """
        first = node == 0
        eps, i1, i2 = self.eps, self.i1, self.i2
        rootI2, xi, intact = self.rootI2, self.xi, self.intact
        alpha, breakage = self.alphaNodal, self.breakageNodal
        trace, floor = self.trace, self.floor
        mu0, lambda0, gammaR, xi0 = self.mu0, self.lambda0, self.gammaR, self.xi0
        aB = self.aB

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
        # solid shear modulus degrades with alpha through gammaR, the granular
        # branch is the aB0..aB3 polynomial in xi.
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
                    3.0 * aB[0]
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
                + aB[3] * yf.mul(xi["l"], i1["l"])
                + aB[1] * rootI2["l"],
            )
            * trace["c"],
        ]

        # The cell as a whole: the means the source guard asks for, and the
        # square of its fastest wave at this node. The square, because that is
        # what the moduli are affine in; the root is taken once per face, by
        # whoever reads the integral.
        speedSquared = self.rhoInv * (lambda0 + twoMuEff["l"])
        statements += [
            self.meanAlpha["u"]
            <= alpha["l"] * self.weights["l"] * self.unitColumn["u"],
            self.meanBreakage["u"]
            <= breakage["l"] * self.weights["l"] * self.unitColumn["u"],
        ]
        nodeSpeed = yf.mul(yf.max(speedSquared, "l"), self.unitColumn["u"])
        statements += [
            (
                self.waveIntegral["u"] <= weight * nodeSpeed
                if first
                else self.waveIntegral["u"]
                <= self.waveIntegral["u"] + weight * nodeSpeed
            ),
            (
                self.intervalLength["u"] <= weight * self.unitColumn["u"]
                if first
                else self.intervalLength["u"]
                <= self.intervalLength["u"] + weight * self.unitColumn["u"]
            ),
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
            self.sourceAlpha["l"]
            <= yf.where(
                growing["l"],
                self.damageRate * drive["l"],
                self.healingRate * drive["l"],
            ),
            self.sourceBreakage["l"]
            <= self.breakageRate
            * yf.mul(
                growing["l"],
                yf.mul(
                    yf.div(
                        1.0,
                        1.0
                        + yf.exp(yf.div(critical["l"] - alpha["l"], self.betaAlpha)),
                    ),
                    drive["l"],
                ),
            ),
        ]

        # What leaves the node: the stress and the source under the quadrature
        # weight, and the two internal variables marched to the next node.
        accSigma, accSource = self.sigmaIntegral, self.sourceIntegral
        alphaColumn = accSource["ln"].subslice("n", 0, 1)
        breakageColumn = accSource["ln"].subslice("n", 1, 2)
        statements += [
            (
                accSigma["lc"] <= weight * sigma["lc"]
                if first
                else accSigma["lc"] <= accSigma["lc"] + weight * sigma["lc"]
            )
        ]
        statements += [
            (
                alphaColumn <= weight * self.sourceAlpha["l"] * self.unitColumn["n"]
                if first
                else accSource["ln"].subslice("n", 0, 1)
                <= accSource["ln"].subslice("n", 0, 1)
                + weight * self.sourceAlpha["l"] * self.unitColumn["n"]
            )
        ]
        statements += [
            (
                breakageColumn
                <= weight * self.sourceBreakage["l"] * self.unitColumn["n"]
                if first
                else accSource["ln"].subslice("n", 1, 2)
                <= accSource["ln"].subslice("n", 1, 2)
                + weight * self.sourceBreakage["l"] * self.unitColumn["n"]
            )
        ]
        statements += [
            alpha["l"] <= alpha["l"] + march * self.sourceAlpha["l"],
            breakage["l"] <= breakage["l"] + march * self.sourceBreakage["l"],
        ]
        return statements

    def finishStatements(self):
        """What the step leaves behind, projected back to modal form once."""
        stress = self.transportGroupSlice("sigma")
        wave = self.transportGroupSlice("waveIntegral")
        interval = self.transportGroupSlice("interval")
        projection = self.db.projectQP[self.t("kl")]
        return [
            self.I["kc"].subslice("c", *stress)
            <= projection * self.sigmaIntegral["lc"],
            self.sourceI["kn"] <= projection * self.sourceIntegral["ln"],
            # A cell value in a modal column is the coefficient of the
            # constant basis function and nothing else.
            self.I["kc"].subslice("c", *wave)
            <= self.constantMode["k"] * self.waveIntegral["c"],
            self.I["kc"].subslice("c", *interval)
            <= self.constantMode["k"] * self.intervalLength["c"],
        ]

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
