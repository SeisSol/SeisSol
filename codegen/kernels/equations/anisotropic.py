# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
# SPDX-FileContributor: Sebastian Wolf

from kernels.equations.elastic import ElasticADERDG
from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile


class AnisotropicADERDG(ElasticADERDG):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir, memLayout)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(f"{matricesDir}/equation-anisotropic.xml", clones)
        )
        memoryLayoutFromFile(memLayout, self.db, clones)

        self.kwargs = kwargs

    def name(self):
        return "anisotropic"


def kernel_class(**kwargs):
    return AnisotropicADERDG
