# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
# @author Jinwen Pan (jinwen.pan AT tum.de)
#

import numpy as np
from yateto.input import parseXMLMatrixFile, memoryLayoutFromFile

from ..aderdg import LinearADERDG


class AcousticADERDG(LinearADERDG):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile("{}/star_acoustic.xml".format(matricesDir), clones)
        )

        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs

    # The 4 quantities are pressure and three velocity components in acoustic materials.
    def numberOfQuantities(self):
        return 4

    def starMatrix(self, dim):
        return self.db.star[dim]

    def addLocal(self, generator, targets):
        super().addLocal(generator, targets)

    def extractVelocities(self):
        extractVelocitiesSPP = np.zeros((3, self.numberOfQuantities()))
        extractVelocitiesSPP[0, 1] = 1
        extractVelocitiesSPP[1, 2] = 1
        extractVelocitiesSPP[2, 3] = 1
        return extractVelocitiesSPP

    def extractTractions(self):
        extractTractionsSPP = np.zeros((3, self.numberOfQuantities()))
        extractTractionsSPP[0, 0] = 1
        return extractTractionsSPP


EQUATION_CLASS = AcousticADERDG
