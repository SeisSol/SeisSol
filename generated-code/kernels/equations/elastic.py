# SPDX-FileCopyrightText: 2016-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.aderdg import LinearADERDG
from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile


class ElasticADERDG(LinearADERDG):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(parseXMLMatrixFile("{}/star.xml".format(matricesDir), clones))

        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs

    def numberOfQuantities(self):
        return 9

    def name(self):
        return "elastic"

    def starMatrix(self, dim):
        return self.db.star[dim]

    def addLocal(self, generator, targets):
        super().addLocal(generator, targets)


EQUATION_CLASS = ElasticADERDG
