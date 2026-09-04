# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from kernels.aderdg.stp import STP
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind
from yateto.input import memoryLayoutFromFile, parseJSONMatrixFile, parseXMLMatrixFile


def choose(n, k):
    num = np.prod(np.arange(n, n - k, -1))
    denom = np.prod(np.arange(1, k + 1))
    return num // denom


class PoroelasticADERDG(STP):
    def __init__(
        self,
        order,
        multipleSimulations,
        matricesDir,
        memLayout,
        numMechanisms,
        **kwargs,
    ):

        super().__init__(
            order, multipleSimulations, matricesDir, memLayout, numMechanisms
        )
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(f"{matricesDir}/equation-poroelastic.xml", clones)
        )
        self.db.update(parseJSONMatrixFile(f"{matricesDir}/stp_{order}.json", clones))

        memoryLayoutFromFile(memLayout, self.db, clones)

        self.kwargs = kwargs

    def primaryGroups(self):
        return [
            QuantityGroup("s", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
            QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
            QuantityGroup("p", QuantityKind.SCALAR, FaceRole.EXTRA_TRACTION),
            QuantityGroup("vf", QuantityKind.VECTOR, FaceRole.EXTRA_VELOCITY),
        ]

    def starMatrix(self, dim):
        return self.db.star[dim]

    def sourceMatrix(self):
        return self.db.ET

    def name(self):
        return "poroelastic"


def kernel_class(**kwargs):
    return PoroelasticADERDG
