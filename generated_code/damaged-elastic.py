#!/usr/bin/env python3

import numpy as np
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile, memoryLayoutFromFile
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

from elastic import ElasticADERDG as ADERDGBase
from multSim import OptionalDimTensor

class DamagedElasticADERDG(ADERDGBase):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
    super().__init__(order, multipleSimulations, matricesDir, memLayout)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update( parseXMLMatrixFile('{}/star.xml'.format(matricesDir), clones) )
    memoryLayoutFromFile(memLayout, self.db, clones)

  def addInit(self, generator):
      super().addInit(generator)

  def add_include_tensors(self, include_tensors):
      super().add_include_tensors(include_tensors)
