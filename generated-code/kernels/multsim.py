# SPDX-FileCopyrightText: 2019-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#

from yateto.type import Tensor
from yateto.ast.node import IndexedTensor
from yateto.memory import DenseMemoryLayout


class OptionalDimTensor(Tensor):
    # dimSize = 1 is considered optional
    def __init__(
        self,
        name,
        optName,
        optSize,
        optPos,
        shape,
        spp=None,
        memoryLayoutClass=DenseMemoryLayout,
        alignStride=False,
    ):
        self._optName = optName
        self._optSize = optSize
        self._optPos = optPos
        shape = self.insertOptDim(shape, (self._optSize,))
        super().__init__(name, shape, spp, memoryLayoutClass, alignStride)

    def hasOptDim(self):
        return self._optSize > 1

    def insertOptDim(self, sliceable, item):
        if self.hasOptDim():
            return sliceable[0 : self._optPos] + item + sliceable[self._optPos :]
        return sliceable

    def __getitem__(self, indexNames):
        indexNames = self.insertOptDim(indexNames, self._optName)
        return IndexedTensor(self, indexNames)

    def optName(self):
        return self._optName

    def optSize(self):
        return self._optSize

    def optPos(self):
        return self._optPos
