# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from kernels.aderdg.stp import STP
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

    def numQuantities(self):
        return 13

    def starMatrix(self, dim):
        return self.db.star[dim]

    def sourceMatrix(self):
        return self.db.ET

    def extractVelocities(self):
        extractVelocitiesSPP = np.zeros((4, self.numQuantities()))
        extractVelocitiesSPP[0, 6] = 1
        extractVelocitiesSPP[1, 7] = 1
        extractVelocitiesSPP[2, 8] = 1
        extractVelocitiesSPP[3, 10] = 1
        return extractVelocitiesSPP

    def extractTractions(self):
        extractTractionsSPP = np.zeros((4, self.numQuantities()))
        extractTractionsSPP[0, 0] = 1
        extractTractionsSPP[1, 3] = 1
        extractTractionsSPP[2, 5] = 1
        extractTractionsSPP[3, 9] = 1
        return extractTractionsSPP

    def name(self):
        return "poroelastic"

    def transformationSpp(self):
        spp = np.zeros((self.numQuantities(), self.numQuantities()), dtype=bool)
        spp[0:6, 0:6] = 1
        spp[6:9, 6:9] = 1
        spp[9, 9] = 1
        spp[10:13, 10:13] = 1
        return spp

    def transformationInvSpp(self):
        return self.transformationSpp()


def kernel_class(**kwargs):
    return PoroelasticADERDG
