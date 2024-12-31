# SPDX-FileCopyrightText: 2016-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
# @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
#

from yateto.input import parseXMLMatrixFile, memoryLayoutFromFile

from kernels.equations.elastic import ElasticADERDG as ADERDGBase


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


EQUATION_CLASS = AnisotropicADERDG
