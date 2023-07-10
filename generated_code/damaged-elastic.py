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
    self.db.update( parseJSONMatrixFile(f'{matricesDir}/dr_stroud_matrices_{self.order}.json', clones, alignStride=self.alignStride, transpose=self.transpose) )

    # For storing nodal values
    qNShape = (numberOfNodes, self.numberOfQuantities())
    self.QNodal = OptionalDimTensor('QNodal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qNShape, alignStride=True)
    dqMShape = (numberOfNodes, self.numberOfQuantities())
    self.dQModal = OptionalDimTensor('dQModal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), dqMShape, alignStride=True)
    FNShape = (numberOfNodes, self.numberOfQuantities())
    self.FNodal = OptionalDimTensor('FNodal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), FNShape, alignStride=True)

    # trans_inv_spp_T = self.transformation_inv_spp().transpose()
    # self.TinvT = Tensor('TinvT', trans_inv_spp_T.shape, spp=trans_inv_spp_T)
    trans_spp_T = self.transformation_spp().transpose()
    self.TT = Tensor('TT', trans_spp_T.shape, spp=trans_spp_T)
    trans_inv_spp_T = self.transformation_inv_spp().transpose()
    self.TinvT = Tensor('TinvT', trans_inv_spp_T.shape, spp=trans_inv_spp_T)


    # WShape = (self.order,)
    # self.Weights = OptionalDimTensor('Weights', 's', multipleSimulations, 0, WShape, alignStride=True)

    memoryLayoutFromFile(memLayout, self.db, clones)

  def numberOfQuantities(self):
    return 10

  def addInit(self, generator):
    super().addInit(generator)

    # For cell average
    generator.add('cellAve', self.QAve['p'] <= self.db.phiAve[self.t('l')] * self.Q['lp'] * 6.0 )
    generator.add('transposeTRot', self.TT['ij'] <= self.T['ji'])


  def addTime(self, generator, targets):
    super().addTime(generator, targets)

    generator.add('damageConvertToNodal', self.QNodal['kp'] <= self.db.v[self.t('kl')] * self.Q['lp'] )
    generator.add('damageAssignFToDQ', self.dQModal['kp'] <= self.db.vInv[self.t('kl')] * self.FNodal['lp'])

    qNShape = (self.t(self.db.v.shape())[0], self.numberOfQuantities())
    Tweight = Scalar('Tweight')
    for i in range(0,self.order):
      QTNodal = OptionalDimTensor('QTNodal({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qNShape, alignStride=True)
      generator.add(f'damageIntegration({i})', self.Q['kp'] <= self.Q['kp'] + self.db.vInv[self.t('kl')] * QTNodal['lp'] * Tweight )

    # Kernels for volumetric flux integration
    gradXiEtaZetaX0 = Scalar('gradXiEtaZetaX0')
    gradXiEtaZetaY0 = Scalar('gradXiEtaZetaY0')
    gradXiEtaZetaZ0 = Scalar('gradXiEtaZetaZ0')
    gradXiEtaZetaX1 = Scalar('gradXiEtaZetaX1')
    gradXiEtaZetaY1 = Scalar('gradXiEtaZetaY1')
    gradXiEtaZetaZ1 = Scalar('gradXiEtaZetaZ1')
    gradXiEtaZetaX2 = Scalar('gradXiEtaZetaX2')
    gradXiEtaZetaY2 = Scalar('gradXiEtaZetaY2')
    gradXiEtaZetaZ2 = Scalar('gradXiEtaZetaZ2')
    # TweightN = Scalar('TweightN')

    numberOfNodes = self.t(self.db.v.shape())[0]
    fluxVolShape = (numberOfNodes, self.numberOfQuantities())

    for i in range(0,self.order):
      FluxVolX = OptionalDimTensor('FluxVolX({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), fluxVolShape, alignStride=True)
      FluxVolY = OptionalDimTensor('FluxVolY({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), fluxVolShape, alignStride=True)
      FluxVolZ = OptionalDimTensor('FluxVolZ({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), fluxVolShape, alignStride=True)

      volumeNonl = (self.Q['kp'] <= \
                    self.Q['kp'] + \
                    ( \
                      self.db.kDivM[0][self.t('kl')] \
                      * self.db.vInv[self.t('lm')] \
                      * (gradXiEtaZetaX0 * FluxVolX['mp']  \
                  + gradXiEtaZetaY0 * FluxVolY['mp'] \
                  + gradXiEtaZetaZ0 * FluxVolZ['mp'] ) \
                    + self.db.kDivM[1][self.t('kl')] \
                      * self.db.vInv[self.t('lm')] \
                      * (gradXiEtaZetaX1 * FluxVolX['mp']  \
                  + gradXiEtaZetaY1 * FluxVolY['mp'] \
                  + gradXiEtaZetaZ1 * FluxVolZ['mp']) \
                    + self.db.kDivM[2][self.t('kl')] \
                      * self.db.vInv[self.t('lm')] \
                      * (gradXiEtaZetaX2 * FluxVolX['mp']  \
                  + gradXiEtaZetaY2 * FluxVolY['mp'] \
                  + gradXiEtaZetaZ2 * FluxVolZ['mp']) \
                    ) #* TweightN
                  )
      generator.add(f'nonlinearVolumeIntegration({i})', volumeNonl)

      # Kernel for Interpolating on the surface quadrature points
      numberOfPoints = self.db.resample.shape()[0]
      gShape = (numberOfPoints, self.numberOfQuantities())
      QInterpolated = OptionalDimTensor('QInterpolated', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), gShape, alignStride=True)

      def interpolateQGenerator(i,h):
        return QInterpolated['kp'] <= self.db.V3mTo2n[i,h][self.t('kl')] * self.Q['lp'] # * self.TinvT['qp']

      # interpolateQPrefetch = lambda i,h: QInterpolated
      generator.addFamily(f'nonlEvaluateAndRotateQAtInterpolationPoints',
                          simpleParameterSpace(4,4),
                          interpolateQGenerator)

      # Surface integration based on the Rusanov flux values on surface quadrature points.
      fluxScale = Scalar('fluxScale')
      Flux = OptionalDimTensor('Flux', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), gShape, alignStride=True)

      nodalFluxGenerator = lambda i,h: self.extendedQTensor()['kp'] <= self.extendedQTensor()['kp'] + fluxScale * self.db.V3mTo2nTWDivM[i,h][self.t('kl')] * Flux['lp'] # * self.TT['qp']
      # nodalFluxPrefetch = lambda i,h: self.I

      generator.addFamily(f'nonlinearSurfaceIntegral',
                          simpleParameterSpace(4,4),
                          nodalFluxGenerator)

  def add_include_tensors(self, include_tensors):
    super().add_include_tensors(include_tensors)
