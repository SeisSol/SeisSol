# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
# SPDX-FileContributor: Sebastian Wolf

import numpy as np
from kernels.equations.elastic import ElasticADERDG as ADERDGBase
from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile


class AnisotropicADERDG(ADERDGBase):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir, memLayout)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile("{}/star_anisotropic.xml".format(matricesDir), clones)
        )
        memoryLayoutFromFile(memLayout, self.db, clones)

        self.kwargs = kwargs

    def name(self):
        return "anisotropic"

    def addInit(self, generator):
        super().addInit(generator)

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)

    def extractTractions(self):
        extractTractionsSPP = np.zeros((3, self.numberOfQuantities()))
        extractTractionsSPP[0, 0] = 1
        extractTractionsSPP[1, 0] = 1
        extractTractionsSPP[2, 0] = 1
        extractTractionsSPP[0, 3] = 1
        extractTractionsSPP[1, 3] = 1
        extractTractionsSPP[2, 3] = 1
        extractTractionsSPP[0, 5] = 1
        extractTractionsSPP[1, 5] = 1
        extractTractionsSPP[2, 5] = 1
        return extractTractionsSPP


EQUATION_CLASS = AnisotropicADERDG
