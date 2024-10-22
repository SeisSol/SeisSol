# SPDX-FileCopyrightText: 2019-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
# @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
#
  
import numpy as np
from yateto import Tensor, Scalar
from kernels.multsim import OptionalDimTensor

def addKernels(generator, aderdg):
  numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
  numberOfQuantities = aderdg.numberOfQuantities()
  order = aderdg.order
  ## Point sources
  mStiffnessTensor = Tensor('stiffnessTensor', (3,3,3,3))
  mSlip = Tensor('mSlip', (3,))
  mNormal = Tensor('mNormal', (3,))
  mArea = Scalar('mArea')
  basisFunctionsAtPoint = Tensor('basisFunctionsAtPoint', (numberOf3DBasisFunctions,))
  basisFunctionDerivativesAtPoint = Tensor('basisFunctionDerivativesAtPoint', (numberOf3DBasisFunctions, 3))
  timeBasisFunctionsAtPoint = Tensor('timeBasisFunctionsAtPoint', (order,))
  mInvJInvPhisAtSources = Tensor('mInvJInvPhisAtSources', (numberOf3DBasisFunctions,))
  JInv = Scalar('JInv')

  generator.add('computeMInvJInvPhisAtSources',
    mInvJInvPhisAtSources['k'] <= JInv * aderdg.db.M3inv['kl'] * basisFunctionsAtPoint['l'])

  #extract the moment tensors entries in SeisSol ordering (xx, yy, zz, xy, yz, xz)
  assert(numberOfQuantities >= 6)
  momentToNRF_spp = np.zeros((numberOfQuantities, 3, 3))
  momentToNRF_spp[0, 0, 0] = 1
  momentToNRF_spp[1, 1, 1] = 1
  momentToNRF_spp[2, 2, 2] = 1
  momentToNRF_spp[3, 0, 1] = 1
  momentToNRF_spp[4, 1, 2] = 1
  momentToNRF_spp[5, 0, 2] = 1
  momentToNRF = Tensor('momentToNRF', (numberOfQuantities, 3, 3), spp=momentToNRF_spp) 

  momentNRFKernel = momentToNRF['tpq'] * mArea * mStiffnessTensor['pqij'] * mSlip['i'] * mNormal['j'] 

  if aderdg.Q.hasOptDim():
    sourceNRF = aderdg.Q['kt'] <= aderdg.Q['kt'] + mInvJInvPhisAtSources['k'] * momentNRFKernel * aderdg.oneSimToMultSim['s'] 
  else:
    sourceNRF = aderdg.Q['kt'] <= aderdg.Q['kt'] + mInvJInvPhisAtSources['k'] * momentNRFKernel 
  generator.add('sourceNRF', sourceNRF)

  momentFSRM = Tensor('momentFSRM', (numberOfQuantities,))
  stfIntegral = Scalar('stfIntegral')
  if aderdg.Q.hasOptDim():
    sourceFSRM = aderdg.Q['kp'] <= aderdg.Q['kp'] + stfIntegral * mInvJInvPhisAtSources['k'] * momentFSRM['p'] * aderdg.oneSimToMultSim['s']
  else:
    sourceFSRM = aderdg.Q['kp'] <= aderdg.Q['kp'] + stfIntegral * mInvJInvPhisAtSources['k'] * momentFSRM['p']
  generator.add('sourceFSRM', sourceFSRM)

  ## Receiver output
  QAtPoint = OptionalDimTensor('QAtPoint', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfQuantities,))
  evaluateDOFSAtPoint = QAtPoint['p'] <= aderdg.Q['kp'] * basisFunctionsAtPoint['k']
  generator.add('evaluateDOFSAtPoint', evaluateDOFSAtPoint)
  QDerivativeAtPoint = OptionalDimTensor('QDerivativeAtPoint', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfQuantities, 3))
  evaluateDerivativeDOFSAtPoint = QDerivativeAtPoint['pd'] <= aderdg.Q['kp'] * basisFunctionDerivativesAtPoint['kd']
  generator.add('evaluateDerivativeDOFSAtPoint', evaluateDerivativeDOFSAtPoint)

  stpShape = (numberOf3DBasisFunctions, numberOfQuantities, order)
  spaceTimePredictor = OptionalDimTensor('spaceTimePredictor', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), stpShape, alignStride=True)
  evaluateDOFSAtPointSTP = QAtPoint['p'] <= spaceTimePredictor['kpt'] * basisFunctionsAtPoint['k'] * timeBasisFunctionsAtPoint['t']
  generator.add('evaluateDOFSAtPointSTP', evaluateDOFSAtPointSTP)
  spaceTimePredictor = OptionalDimTensor('spaceTimePredictor', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), stpShape, alignStride=True)
  evaluateDerivativeDOFSAtPointSTP = QDerivativeAtPoint['pd'] <= spaceTimePredictor['kpt'] * basisFunctionDerivativesAtPoint['kd'] * timeBasisFunctionsAtPoint['t']
  generator.add('evaluateDerivativeDOFSAtPointSTP', evaluateDerivativeDOFSAtPointSTP)
