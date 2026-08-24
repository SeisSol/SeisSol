# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile

from .viscoelastic2 import Viscoelastic2ADERDG


class Viscoacoustic2ADERDG(Viscoelastic2ADERDG):
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
            order,
            multipleSimulations,
            matricesDir,
            memLayout,
            numMechanisms,
            **kwargs,
        )

        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(
                "{}/matrices_viscoacoustic.xml".format(matricesDir), clones
            )
        )

        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs

    def numQuantities(self):
        return 4

    def numAnelasticQuantities(self):
        return 1


EQUATION_CLASS = Viscoacoustic2ADERDG
