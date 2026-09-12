# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Continuum damage-breakage rheology.

The state is a strain tensor, a velocity, and the two internal variables of the
rheology: the damage alpha and the breakage B. Neither internal variable takes
part in the face coupling -- they evolve through source terms that are local to
the cell -- so both sit in the layout without a face role.

The strain carries the traction role even though it is a strain rather than a
stress. The two transform alike under the face rotation, which is what the role
drives today; that the mechanical traction is a derived quantity here, obtained
from sigma(eps, alpha, B), is a distinction the face machinery does not yet
draw.
"""

from kernels.aderdg.nonlinearck import NonLinearCK
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind


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


def kernel_class(**kwargs):
    solver = kwargs["solver"].lower()
    if solver == "nonlinearck":
        return DamageADERDG
    raise NotImplementedError(f"{solver} cannot advance a material with damage.")
