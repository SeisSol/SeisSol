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
                growing["l"] <= yf.greater(xi["l"] + xi0, floor),
                sourceAlpha["l"]
                <= yf.where(
                    growing["l"], damageRate * drive["l"], healingRate * drive["l"]
                ),
                sourceBreakage["l"]
                <= breakageRate * yf.mul(growing["l"], yf.mul(switch["l"], drive["l"])),
            ],
        )


def kernel_class(**kwargs):
    solver = kwargs["solver"].lower()
    if solver == "nonlinearck":
        return DamageADERDG
    raise NotImplementedError(f"{solver} cannot advance a material with damage.")
