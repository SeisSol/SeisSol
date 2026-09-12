# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Jinwen Pan

from kernels.aderdg.linearck import LinearCK
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind


class AcousticADERDG(LinearCK):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir)
        self.configure(matricesDir, memLayout, kwargs)

    # The 4 quantities are pressure and three velocity components
    # in acoustic materials.
    def primaryGroups(self):
        return [
            QuantityGroup("pprime", QuantityKind.SCALAR, FaceRole.TRACTION),
            QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
        ]

    def name(self):
        return "acoustic"


def kernel_class(**kwargs):
    return AcousticADERDG
