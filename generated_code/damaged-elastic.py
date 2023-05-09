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
    # self.db.update( parseXMLMatrixFile('{}/star.xml'.format(matricesDir), clones) )
    self.db.update( parseXMLMatrixFile('{}/matrices_damaged-elastic.xml'.format(matricesDir), clones) )
    self.db.update( parseXMLMatrixFile('{}/plasticity_{}_matrices_{}.xml'.format(matricesDir, 'nb', self.order), clones=dict(), alignStride=self.alignStride) )
    numberOfNodes = self.t(self.db.v.shape())[0]
    # For cell-average
    self.db.update(parseXMLMatrixFile('{}/phi_ave_{}.xml'.format(matricesDir, order), transpose=self.transpose, alignStride=self.alignStride))
    
    qAveShape = (self.numberOfQuantities(),)
    self.QAve = OptionalDimTensor('QAve', 's', multipleSimulations, 0, qAveShape, alignStride=True)

    # For storing nodal values
    qNShape = (numberOfNodes, self.numberOfQuantities())
    self.QNodal = OptionalDimTensor('QNodal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qNShape, alignStride=True)
    dqMShape = (numberOfNodes, self.numberOfQuantities())
    self.dQModal = OptionalDimTensor('dQModal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), dqMShape, alignStride=True)
    FNShape = (numberOfNodes, self.numberOfQuantities())
    self.FNodal = OptionalDimTensor('FNodal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), FNShape, alignStride=True)
    # WShape = (self.order,)
    # self.Weights = OptionalDimTensor('Weights', 's', multipleSimulations, 0, WShape, alignStride=True)

    memoryLayoutFromFile(memLayout, self.db, clones)

  def numberOfQuantities(self):
    return 10

  def addInit(self, generator):
    super().addInit(generator)

    # For cell average
    generator.add('cellAve', self.QAve['p'] <= self.db.phiAve[self.t('l')] * self.Q['lp'] * 6.0 )


  def addTime(self, generator, targets):
    super().addTime(generator, targets)

    generator.add('damageConvertToNodal', self.QNodal['kp'] <= self.db.v[self.t('kl')] * self.Q['lp'] )
    generator.add('damageAssignFToDQ', self.dQModal['kp'] <= self.db.vInv[self.t('kl')] * self.FNodal['lp'])
    
    qNShape = (self.t(self.db.v.shape())[0], self.numberOfQuantities())
    Tweight = Scalar('Tweight')
    for i in range(0,self.order):
      QTNodal = OptionalDimTensor('QTNodal({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qNShape, alignStride=True)
      generator.add(f'damageIntegration({i})', self.Q['kp'] <= self.Q['kp'] + self.db.vInv[self.t('kl')] * QTNodal['lp'] * Tweight )
    

  def add_include_tensors(self, include_tensors):
    super().add_include_tensors(include_tensors)
