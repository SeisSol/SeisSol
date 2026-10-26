# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.aderdg.linearck import LinearCK
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind


class ElasticADERDG(LinearCK):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir)
        self.configure(matricesDir, memLayout, kwargs)

    def primaryGroups(self):
        return [
            QuantityGroup("s", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
            QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
        ]

    def name(self):
        return "elastic"


def kernel_class(**kwargs):
    return ElasticADERDG
