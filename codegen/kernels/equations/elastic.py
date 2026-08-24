# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels.aderdg.linearck import LinearCK
from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile


class ElasticADERDG(LinearCK):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(f"{matricesDir}/equation-elastic.xml", clones)
        )

        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs

    def numQuantities(self):
        return 9

    def name(self):
        return "elastic"

    def starMatrix(self, dim):
        return self.db.star[dim]

    def extractVelocities(self):
        extractVelocitiesSPP = np.zeros((3, self.numQuantities()))
        extractVelocitiesSPP[0, 6] = 1
        extractVelocitiesSPP[1, 7] = 1
        extractVelocitiesSPP[2, 8] = 1
        return extractVelocitiesSPP

    def extractTractions(self):
        extractTractionsSPP = np.zeros((3, self.numQuantities()))
        extractTractionsSPP[0, 0] = 1
        extractTractionsSPP[1, 3] = 1
        extractTractionsSPP[2, 5] = 1
        return extractTractionsSPP


def kernel_class(**kwargs):
    return ElasticADERDG
