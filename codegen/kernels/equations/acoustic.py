# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Jinwen Pan

import numpy as np
from kernels.aderdg.linearck import LinearCK
from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind


class AcousticADERDG(LinearCK):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(f"{matricesDir}/equation-acoustic.xml", clones)
        )

        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs

    # The 4 quantities are pressure and three velocity components
    # in acoustic materials.
    def primaryGroups(self):
        return [
            QuantityGroup("pprime", QuantityKind.SCALAR, FaceRole.TRACTION),
            QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
        ]

    def starMatrix(self, dim):
        return self.db.star[dim]

    def name(self):
        return "acoustic"


def kernel_class(**kwargs):
    return AcousticADERDG
