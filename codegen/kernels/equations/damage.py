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

    def primaryGroups(self):
        return [
            QuantityGroup("eps", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
            QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
            QuantityGroup("alpha", QuantityKind.SCALAR),
            QuantityGroup("breakage", QuantityKind.SCALAR),
        ]

    def name(self):
        return "damage"

    def addConstitutive(self, generator):
        nodes = self.num3DQuadraturePoints()
        nq = self.numQuantities()

        lambda0 = Scalar("lambda0")
        mu0 = Scalar("mu0")
        gammaR = Scalar("gammaR")
        xi0 = Scalar("xi0")
        damageRate = Scalar("damageRate")
        breakageRate = Scalar("breakageRate")
        healingRate = Scalar("healingRate")
        betaAlpha = Scalar("betaAlpha")
        aB = [Scalar(f"aB{i}") for i in range(4)]

        # Guards xi against a vanishing second invariant. Tied to the working
        # precision rather than to the model, hence not a model parameter.
        floor = Scalar("invariantFloor")

        ones = Tensor("ones", (nodes,), np.ones(nodes))
        trace = Tensor("traceSelect", (6,), TRACE)
        voigt = Tensor("voigtSquare", (6,), VOIGT_SQUARE)
        pickAlpha = Tensor("selectAlpha", (nq,), unit(ALPHA, nq))
        pickBreakage = Tensor("selectBreakage", (nq,), unit(BREAKAGE, nq))
        epsInit = Tensor("epsInit", (6,))
        unitColumn = Tensor("unitColumn", (1,), np.ones(1))

        eps = self.nodalTensor("epsTotal", 6)
        epsSquare = self.nodalTensor("epsSquare", 6)
        i1 = self.nodalTensor("invariantI1")
        i2 = self.nodalTensor("invariantI2")
        rootI2 = self.nodalTensor("rootI2")
        xi = self.nodalTensor("xi")
        alpha = self.nodalTensor("alphaNodal")
        breakage = self.nodalTensor("breakageNodal")
        intact = self.nodalTensor("intact")

        generator.add(
            "damageInvariants",
            [
                eps["lc"]
                <= yf.add(self.QNodal["lc"].subslice("c", 0, 6), epsInit["c"]),
                epsSquare["lc"] <= yf.mul(eps["lc"], eps["lc"]),
                i1["l"] <= eps["lc"] * trace["c"],
                i2["l"] <= epsSquare["lc"] * voigt["c"],
                rootI2["l"] <= yf.sqrt(yf.maximum(i2["l"], floor)),
                xi["l"]
                <= yf.where(
                    yf.greater(i2["l"], floor),
                    yf.div(i1["l"], rootI2["l"]),
                    0.0,
                ),
                alpha["l"] <= self.QNodal["lp"] * pickAlpha["p"],
                breakage["l"] <= self.QNodal["lp"] * pickBreakage["p"],
                intact["l"] <= 1.0 - breakage["l"],
            ],
        )

        twoMuEff = self.nodalTensor("twoMuEff")
        isoSolid = self.nodalTensor("isotropicSolid")
        isoGranular = self.nodalTensor("isotropicGranular")
        shearGranular = self.nodalTensor("shearGranular")
        sigmaSolid = self.nodalTensor("sigmaSolid", 6)
        sigmaGranular = self.nodalTensor("sigmaGranular", 6)
        sigma = self.nodalTensor("sigmaNodal", 6)

        generator.add(
            "damageStress",
            [
                twoMuEff["l"]
                <= 2.0 * mu0
                - 2.0 * gammaR * xi0 * alpha["l"]
                - gammaR * yf.mul(alpha["l"], xi["l"]),
                isoSolid["l"]
                <= lambda0 * i1["l"] - gammaR * yf.mul(alpha["l"], rootI2["l"]),
                sigmaSolid["lc"] <= yf.mul(twoMuEff["l"], eps["lc"]),
                sigmaSolid["lc"] <= sigmaSolid["lc"] + isoSolid["l"] * trace["c"],
                isoGranular["l"]
                <= 2.0 * aB[2] * i1["l"]
                + aB[3] * yf.mul(xi["l"], i1["l"])
                + aB[1] * rootI2["l"],
                shearGranular["l"]
                <= 3.0 * aB[0]
                + aB[1] * xi["l"]
                - aB[3] * yf.mul(xi["l"], yf.mul(xi["l"], xi["l"])),
                sigmaGranular["lc"] <= yf.mul(shearGranular["l"], eps["lc"]),
                sigmaGranular["lc"]
                <= sigmaGranular["lc"] + isoGranular["l"] * trace["c"],
                sigma["lc"]
                <= yf.add(
                    yf.mul(intact["l"], sigmaSolid["lc"]),
                    yf.mul(breakage["l"], sigmaGranular["lc"]),
                ),
            ],
        )
        self.sigmaNodal = sigma

        rhoInv = Scalar("rhoInv")
        velocity = self.nodalTensor("velocityNodal", 3)
        flux = [self.nodalTensor(f"fluxNodal{axis}", nq) for axis in "XYZ"]
        toFluxV = [
            Tensor(f"velocityToFlux{axis}", (3, nq), fluxMap(VELOCITY_FLUX[d], 3, nq))
            for d, axis in enumerate("XYZ")
        ]
        toFluxS = [
            Tensor(f"stressToFlux{axis}", (6, nq), fluxMap(STRESS_FLUX[d], 6, nq))
            for d, axis in enumerate("XYZ")
        ]

        assembly = [velocity["lm"] <= self.QNodal["lm"].subslice("m", 6, 9)]
        for d in range(3):
            assembly.append(
                flux[d]["lp"]
                <= velocity["lm"] * toFluxV[d]["mp"]
                + rhoInv * self.sigmaNodal["lc"] * toFluxS[d]["cp"]
            )
        generator.add("damageFlux", assembly)
        self.fluxNodal = flux

        # Stage D: the two reductions the cell needs as a whole.
        #
        # The source guard asks whether the cell still has room to damage,
        # which is a property of the cell and not of a node, so alpha and B
        # enter it through their means. The Rusanov dissipation needs one wave
        # speed per cell, and the largest one over the nodes is the safe pick.
        weights = Tensor("quadratureWeights", (nodes,))
        meanAlpha = Tensor("meanAlpha", (1,))
        meanBreakage = Tensor("meanBreakage", (1,))
        waveSpeed = self.nodalTensor("waveSpeedNodal")
        maxWaveSpeed = Tensor("maxWaveSpeed", (1,))

        generator.add(
            "damageCellState",
            [
                meanAlpha["u"] <= alpha["l"] * weights["l"] * unitColumn["u"],
                meanBreakage["u"] <= breakage["l"] * weights["l"] * unitColumn["u"],
                waveSpeed["l"] <= yf.sqrt(rhoInv * (lambda0 + twoMuEff["l"])),
                maxWaveSpeed["u"]
                <= yf.mul(yf.max(waveSpeed["l"], "l"), unitColumn["u"]),
            ],
        )

        quadA = self.nodalTensor("criticalA")
        quadB = self.nodalTensor("criticalB")
        quadC = self.nodalTensor("criticalC")
        rootTerm = self.nodalTensor("criticalRoot")
        critical1 = self.nodalTensor("criticalFromRoot")
        critical2 = self.nodalTensor("criticalFromModuli")
        critical = self.nodalTensor("criticalDamage")
        switch = self.nodalTensor("granularSwitch")
        drive = self.nodalTensor("damageDrive")
        growing = self.nodalTensor("damageGrowing")
        sourceAlpha = self.nodalTensor("sourceAlpha")
        sourceBreakage = self.nodalTensor("sourceBreakage")

        generator.add(
            "damageSource",
            [
                quadA["l"]
                <= 3.0 * gammaR * gammaR * yf.mul(xi["l"], xi["l"])
                - 3.0 * gammaR * gammaR
                + 6.0 * gammaR * gammaR * xi0 * xi["l"]
                + 4.0 * gammaR * gammaR * xi0 * xi0,
                quadB["l"]
                <= -(8.0 * mu0 + 6.0 * lambda0) * gammaR * xi0
                - gammaR * lambda0 * yf.mul(xi["l"], yf.mul(xi["l"], xi["l"]))
                - 6.0 * gammaR * mu0 * xi["l"],
                quadC["l"] <= (4.0 * mu0 * mu0 + 6.0 * mu0 * lambda0) * ones["l"],
                rootTerm["l"]
                <= yf.sqrt(
                    yf.maximum(
                        yf.mul(quadB["l"], quadB["l"])
                        - 4.0 * yf.mul(quadA["l"], quadC["l"]),
                        floor,
                    )
                ),
                critical1["l"] <= yf.div(-quadB["l"] - rootTerm["l"], 2.0 * quadA["l"]),
                critical2["l"] <= yf.div(2.0 * mu0, gammaR * (xi["l"] + 2.0 * xi0)),
                critical["l"]
                <= yf.minimum(
                    yf.minimum(
                        yf.where(
                            yf.greater(critical1["l"], floor),
                            critical1["l"],
                            1.0,
                        ),
                        yf.where(
                            yf.greater(critical2["l"], floor),
                            critical2["l"],
                            1.0,
                        ),
                    ),
                    1.0,
                ),
                switch["l"]
                <= yf.div(
                    1.0,
                    1.0 + yf.exp(yf.div(critical["l"] - alpha["l"], betaAlpha)),
                ),
                drive["l"]
                <= gammaR * yf.mul(intact["l"], yf.mul(i2["l"], xi["l"] + xi0)),
                growing["l"]
                <= yf.logical_and(
                    yf.greater(xi["l"] + xi0, floor),
                    yf.logical_and(
                        yf.less(yf.sum(meanAlpha["u"], "u"), 1.0),
                        yf.less(yf.sum(meanBreakage["u"], "u"), 1.0),
                    ),
                ),
                sourceAlpha["l"]
                <= yf.where(
                    growing["l"], damageRate * drive["l"], healingRate * drive["l"]
                ),
                sourceBreakage["l"]
                <= breakageRate * yf.mul(growing["l"], yf.mul(switch["l"], drive["l"])),
            ],
        )

    def addFaceFlux(self, generator):
        """The Rusanov flux at the face nodes.

        Both sides are read from their own transported stress rather than
        rebuilt from the other cell's state, so neither side ever needs the
        other's material. The only thing crossing the face besides the state
        and the stress is one wave speed, and taking the larger of the two is
        what keeps the scheme stable.

        The directional maps from the flux assembly are reused, contracted with
        the face normal instead of applied one direction at a time. They carry
        the sign of the flux, so the combination below reads as the plain
        average it is.
        """
        nq = self.numQuantities()
        rhoInv = Scalar("rhoInv")
        rhoInvNeighbor = Scalar("rhoInvNeighbor")
        lambdaMax = Scalar("lambdaMax")

        normal = Tensor("faceNormal", (3,))
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

        sigmaFace = self.faceTensor("sigmaAtFace", 6)
        sigmaFaceNeighbor = self.faceTensor("sigmaAtFaceNeighbor", 6)
        fluxLocal = self.faceTensor("fluxAtFaceLocal", nq)
        fluxNeighbor = self.faceTensor("fluxAtFaceNeighbor", nq)

        generator.add(
            "damageRusanov",
            [
                fluxLocal["kp"]
                <= self.QAtFace["km"].subslice("m", 6, 9)
                * velocityMap["dmp"]
                * normal["d"]
                + rhoInv * sigmaFace["kc"] * stressMap["dcp"] * normal["d"],
                fluxNeighbor["kp"]
                <= self.QAtFaceNeighbor["km"].subslice("m", 6, 9)
                * velocityMap["dmp"]
                * normal["d"]
                + rhoInvNeighbor
                * sigmaFaceNeighbor["kc"]
                * stressMap["dcp"]
                * normal["d"],
                self.fluxAtFace["kp"]
                <= 0.5 * (fluxLocal["kp"] + fluxNeighbor["kp"])
                - 0.5 * lambdaMax * (self.QAtFaceNeighbor["kp"] - self.QAtFace["kp"]),
            ],
        )


def kernel_class(**kwargs):
    solver = kwargs["solver"].lower()
    if solver == "nonlinearck":
        return DamageADERDG
    raise NotImplementedError(f"{solver} cannot advance a material with damage.")
