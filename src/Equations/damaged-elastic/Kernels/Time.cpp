/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Time kernel of SeisSol.
 **/

#include "Kernels/TimeBase.h"
#include "Kernels/Time.h"
#include "Kernels/GravitationalFreeSurfaceBC.h"
#include <DynamicRupture/Misc.h>
#include <DynamicRupture/Typedefs.hpp>
#include <Initializer/Parameters/ModelParameters.h>
#include <Kernels/precision.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>

#include <cstring>
#include <cassert>
#include <stdint.h>
#include <omp.h>

#include <yateto.h>

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

seissol::kernels::TimeBase::TimeBase() {
  m_derivativesOffsets[0] = 0;
  for (int order = 0; order < CONVERGENCE_ORDER; ++order) {
    if (order > 0) {
      m_derivativesOffsets[order] = tensor::dQ::size(order - 1) + m_derivativesOffsets[order - 1];
    }
  }
}

void seissol::kernels::TimeBase::checkGlobalData(GlobalData const* global, size_t alignment) {
  assert(((uintptr_t)global->stiffnessMatricesTransposed(0)) % alignment == 0);
  assert(((uintptr_t)global->stiffnessMatricesTransposed(1)) % alignment == 0);
  assert(((uintptr_t)global->stiffnessMatricesTransposed(2)) % alignment == 0);
}

void seissol::kernels::Time::setHostGlobalData(GlobalData const* global) {
#ifdef USE_STP
  // Note: We could use the space time predictor for elasticity.
  // This is not tested and experimental
  for (int n = 0; n < CONVERGENCE_ORDER; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        m_krnlPrototype.kDivMTSub(d, n) = init::kDivMTSub::Values[tensor::kDivMTSub::index(d, n)];
      }
    }
    m_krnlPrototype.selectModes(n) = init::selectModes::Values[tensor::selectModes::index(n)];
  }
  m_krnlPrototype.Zinv = init::Zinv::Values;
  m_krnlPrototype.timeInt = init::timeInt::Values;
  m_krnlPrototype.wHat = init::wHat::Values;
#else // USE_STP
  checkGlobalData(global, ALIGNMENT);

  m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;

  projectDerivativeToNodalBoundaryRotated.V3mTo2nFace = global->V3mTo2nFace;

#endif // USE_STP
}

void seissol::kernels::Time::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  checkGlobalData(global.onDevice, deviceAlignment);
  deviceKrnlPrototype.kDivMT = global.onDevice->stiffnessMatricesTransposed;
  deviceDerivativeToNodalBoundaryRotated.V3mTo2nFace = global.onDevice->V3mTo2nFace;

#endif
}

void seissol::kernels::Time::computeAder(double i_timeStepWidth,
                                         LocalData& data,
                                         LocalTmp& tmp,
                                         real o_timeIntegrated[tensor::I::size()],
                                         real* o_timeDerivatives,
                                         bool updateDisplacement) {

  assert(reinterpret_cast<uintptr_t>(data.dofs) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(o_timeIntegrated) % ALIGNMENT == 0);
  assert(o_timeDerivatives == nullptr ||
         reinterpret_cast<uintptr_t>(o_timeDerivatives) % ALIGNMENT == 0);

  // Only a small fraction of cells has the gravitational free surface boundary condition
  updateDisplacement &=
      std::any_of(std::begin(data.cellInformation.faceTypes),
                  std::end(data.cellInformation.faceTypes),
                  [](const FaceType f) { return f == FaceType::freeSurfaceGravity; });

#ifdef USE_STP
  // Note: We could use the space time predictor for elasticity.
  // This is not tested and experimental
  alignas(PAGESIZE_STACK) real stpRhs[tensor::spaceTimePredictor::size()];
  alignas(PAGESIZE_STACK) real stp[tensor::spaceTimePredictor::size()]{};
  kernel::spaceTimePredictor krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }
  krnl.Q = const_cast<real*>(data.dofs);
  krnl.I = o_timeIntegrated;
  krnl.timestep = i_timeStepWidth;
  krnl.spaceTimePredictor = stp;
  krnl.spaceTimePredictorRhs = stpRhs;
  krnl.execute();
#else // USE_STP

#ifdef USE_DAMAGEDELASTIC
  const real damageParameter = data.material.local.Cd;
  const real scalingValue = m_damagedElasticParameters->scalingValue;
  const real breakCoefficient = scalingValue * damageParameter;
  const real betaAlpha = m_damagedElasticParameters->betaAlpha;

  kernel::damageConvertToNodal d_converToKrnl;
  // Compute the nodal solutions
  alignas(PAGESIZE_STACK) real solNData[tensor::QNodal::size()];
  d_converToKrnl.v = init::v::Values;
  d_converToKrnl.QNodal = solNData;
  d_converToKrnl.Q = data.dofs;
  d_converToKrnl.execute();

  // Compute rhs of damage evolution
  alignas(PAGESIZE_STACK) real fNodalData[tensor::FNodal::size()] = {0};
  real* exxNodal = (solNData + 0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* eyyNodal = (solNData + 1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* ezzNodal = (solNData + 2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* exyNodal = (solNData + 3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* eyzNodal = (solNData + 4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* ezxNodal = (solNData + 5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* alphaNodal = (solNData + 9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* breakNodal = (solNData + 10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);

  real alpha_ave = alphaNodal[0];
  real break_ave = breakNodal[0];
  for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS - 1; ++q) {
    break_ave = std::max(break_ave, breakNodal[q]);
    alpha_ave = std::max(alpha_ave, alphaNodal[q]);
  }

  for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q) {
    real EspI, EspII, xi;
    calculateEps(exxNodal, eyyNodal, ezzNodal, exyNodal, eyzNodal, ezxNodal, q, *m_damagedElasticParameters, EspI, EspII, xi);

    // Compute alpha_{cr}
    real aCR = (3.0 * xi * xi - 3.0) * data.material.local.gammaR * data.material.local.gammaR +
               6.0 * xi * data.material.local.gammaR * data.material.local.xi0 *
                   data.material.local.gammaR +
               4.0 * data.material.local.xi0 * data.material.local.gammaR *
                   data.material.local.xi0 * data.material.local.gammaR;

    real bCR = -(8.0 * data.material.local.mu0 + 6.0 * data.material.local.lambda0) *
                   data.material.local.xi0 * data.material.local.gammaR -
               xi * (xi * xi * data.material.local.lambda0 + 6.0 * data.material.local.mu0) *
                   data.material.local.gammaR;

    real cCR = 4.0 * data.material.local.mu0 * data.material.local.mu0 +
               6.0 * data.material.local.mu0 * data.material.local.lambda0;

    real alphaCR1q = (-bCR - std::sqrt(bCR * bCR - 4.0 * aCR * cCR)) / (2.0 * aCR);
    real alphaCR2q = 2.0 * data.material.local.mu0 / data.material.local.gammaR /
                     (xi + 2.0 * data.material.local.xi0);

    real alphaCRq = 1.0;
    if (alphaCR1q > 0.0) {
      if (alphaCR2q > 0.0) {
        alphaCRq = std::min(static_cast<real>(1.0), std::min(alphaCR1q, alphaCR2q));
      }
    }

    if (xi + data.material.local.xi0 > 0) {
      if (alpha_ave < 0.9) {
        if (break_ave < 0.85) {
          fNodalData[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
              (1 - breakNodal[q]) * 1.0 / (std::exp((alphaCRq - alphaNodal[q]) / betaAlpha) + 1.0) *
              breakCoefficient * data.material.local.gammaR * EspII *
              (xi + data.material.local.xi0);
          fNodalData[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
              (1 - breakNodal[q]) * damageParameter * data.material.local.gammaR * EspII *
              (xi + data.material.local.xi0);
        } else {
          fNodalData[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
          fNodalData[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
        }
      } else {
        fNodalData[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
        fNodalData[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
      }
    } else if (alpha_ave > 5e-1) {
      fNodalData[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0 * damageParameter *
                                                              data.material.local.gammaR * EspII *
                                                              (xi + data.material.local.xi0);
      fNodalData[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0 * damageParameter *
                                                               data.material.local.gammaR * EspII *
                                                               (xi + data.material.local.xi0);
    } else {
      fNodalData[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
      fNodalData[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
    }
  }

  // Convert them back to modal space
  alignas(PAGESIZE_STACK) real dQModalData[tensor::dQModal::size()];
  kernel::damageAssignFToDQ d_assignFToDQ;
  d_assignFToDQ.vInv = init::vInv::Values;
  d_assignFToDQ.dQModal = dQModalData;
  d_assignFToDQ.FNodal = fNodalData;
  d_assignFToDQ.execute();

#endif

  alignas(PAGESIZE_STACK) real temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()];
  auto* derivativesBuffer = (o_timeDerivatives != nullptr) ? o_timeDerivatives : temporaryBuffer;

  kernel::derivative krnl = m_krnlPrototype;
  krnl.dQModal = dQModalData;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }

  // Optional source term
  set_ET(krnl, get_ptr_sourceMatrix(data.localIntegration.specific));

  krnl.dQ(0) = const_cast<real*>(data.dofs);
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  intKrnl.dQ(0) = data.dofs;
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }
  // powers in the taylor-series expansion
  intKrnl.power = i_timeStepWidth;
  intKrnl.execute0();

  if (updateDisplacement) {
    // First derivative if needed later in kernel
    std::copy_n(data.dofs, tensor::dQ::size(0), derivativesBuffer);
  } else if (o_timeDerivatives != nullptr) {
    // First derivative is not needed here but later
    // Hence stream it out
    streamstore(tensor::dQ::size(0), data.dofs, derivativesBuffer);
  }

  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
    krnl.execute(der);

    // update scalar for this derivative
    intKrnl.power *= i_timeStepWidth / real(der + 1);
    intKrnl.execute(der);
  }

  // Do not compute it like this if at interface
  // Compute integrated displacement over time step if needed.
  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      if (data.faceDisplacements[face] != nullptr &&
          data.cellInformation.faceTypes[face] == FaceType::freeSurfaceGravity) {
        bc.evaluate(face,
                    projectDerivativeToNodalBoundaryRotated,
                    data.boundaryMapping[face],
                    data.faceDisplacements[face],
                    tmp.nodalAvgDisplacements[face].data(),
                    *this,
                    derivativesBuffer,
                    i_timeStepWidth,
                    data.material,
                    data.cellInformation.faceTypes[face]);
      }
    }
  }
#endif // USE_STP
}

void seissol::kernels::Time::computeBatchedAder(double i_timeStepWidth,
                                                LocalTmp& tmp,
                                                ConditionalPointersToRealsTable& dataTable,
                                                ConditionalMaterialTable& materialTable,
                                                bool updateDisplacement) {
#ifdef ACL_DEVICE
  kernel::gpu_derivative derivativesKrnl = deviceKrnlPrototype;
  kernel::gpu_derivativeTaylorExpansion intKrnl;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    derivativesKrnl.numElements = numElements;
    intKrnl.numElements = numElements;

    intKrnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();

    unsigned starOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      derivativesKrnl.star(i) =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
      derivativesKrnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }

    unsigned derivativesOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derivativesKrnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      intKrnl.dQ(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr());

      derivativesKrnl.extraOffset_dQ(i) = derivativesOffset;
      intKrnl.extraOffset_dQ(i) = derivativesOffset;

      derivativesOffset += tensor::dQ::size(i);
    }

    // stream dofs to the zero derivative
    device.algorithms.streamBatchedData(
        (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
        tensor::Q::Size,
        derivativesKrnl.numElements,
        device.api->getDefaultStream());

    const auto maxTmpMem = yateto::getMaxTmpMemRequired(intKrnl, derivativesKrnl);
    real* tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(maxTmpMem * numElements));

    intKrnl.power = i_timeStepWidth;
    intKrnl.linearAllocator.initialize(tmpMem);
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute0();

    for (unsigned Der = 1; Der < CONVERGENCE_ORDER; ++Der) {
      derivativesKrnl.linearAllocator.initialize(tmpMem);
      derivativesKrnl.streamPtr = device.api->getDefaultStream();
      derivativesKrnl.execute(Der);

      // update scalar for this derivative
      intKrnl.power *= i_timeStepWidth / real(Der + 1);
      intKrnl.linearAllocator.initialize(tmpMem);
      intKrnl.streamPtr = device.api->getDefaultStream();
      intKrnl.execute(Der);
    }
    device.api->popStackMemory();
  }

  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      bc.evaluateOnDevice(face,
                          deviceDerivativeToNodalBoundaryRotated,
                          *this,
                          dataTable,
                          materialTable,
                          i_timeStepWidth,
                          device);
    }
  }
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::flopsAder(unsigned int& o_nonZeroFlops,
                                       unsigned int& o_hardwareFlops) {
  // reset flops
  o_nonZeroFlops = 0;
  o_hardwareFlops = 0;

  // initialization
  o_nonZeroFlops += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for (unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++) {
    o_nonZeroFlops += kernel::derivative::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivative::hardwareFlops(l_derivative);

    // update of time integrated DOFs
    o_nonZeroFlops += kernel::derivativeTaylorExpansion::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(l_derivative);
  }
}

unsigned seissol::kernels::Time::bytesAder() {
  unsigned reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();

  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral(double i_expansionPoint,
                                             double i_integrationStart,
                                             double i_integrationEnd,
                                             const real* i_timeDerivatives,
                                             real o_timeIntegrated[tensor::I::size()]) {
  /*
   * assert alignments.
   */
  assert(((uintptr_t)i_timeDerivatives) % ALIGNMENT == 0);
  assert(((uintptr_t)o_timeIntegrated) % ALIGNMENT == 0);

  // assert that this is a forwared integration in time
  assert(i_integrationStart + (real)1.E-10 > i_expansionPoint);
  assert(i_integrationEnd > i_integrationStart);

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm = (real)1;
  real l_secondTerm = (real)1;
  real l_factorial = (real)1;
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives + m_derivativesOffsets[i];
  }

  // iterate over time derivatives
  for (int der = 0; der < CONVERGENCE_ORDER; ++der) {
    l_firstTerm *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial *= (real)(der + 1);

    intKrnl.power = l_firstTerm - l_secondTerm;
    intKrnl.power /= l_factorial;

    intKrnl.execute(der);
  }
}

void seissol::kernels::Time::computeBatchedIntegral(double i_expansionPoint,
                                                    double i_integrationStart,
                                                    double i_integrationEnd,
                                                    const real** i_timeDerivatives,
                                                    real** o_timeIntegratedDofs,
                                                    unsigned numElements) {
#ifdef ACL_DEVICE
  // assert that this is a forwared integration in time
  assert(i_integrationStart + (real)1.E-10 > i_expansionPoint);
  assert(i_integrationEnd > i_integrationStart);

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real deltaTLower = i_integrationStart - i_expansionPoint;
  real deltaTUpper = i_integrationEnd - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real firstTerm = static_cast<real>(1.0);
  real secondTerm = static_cast<real>(1.0);
  real factorial = static_cast<real>(1.0);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  real* tmpMem = reinterpret_cast<real*>(
      device.api->getStackMemory(intKrnl.TmpMaxMemRequiredInBytes * numElements));

  intKrnl.I = o_timeIntegratedDofs;

  unsigned derivativesOffset = 0;
  for (size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives;
    intKrnl.extraOffset_dQ(i) = derivativesOffset;
    derivativesOffset += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  for (int der = 0; der < CONVERGENCE_ORDER; ++der) {
    firstTerm *= deltaTUpper;
    secondTerm *= deltaTLower;
    factorial *= static_cast<real>(der + 1);

    intKrnl.power = firstTerm - secondTerm;
    intKrnl.power /= factorial;
    intKrnl.linearAllocator.initialize(tmpMem);
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute(der);
  }
  device.api->popStackMemory();
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::computeTaylorExpansion(real time,
                                                    real expansionPoint,
                                                    real const* timeDerivatives,
                                                    real timeEvaluated[tensor::Q::size()]) {
  /*
   * assert alignments.
   */
  assert(((uintptr_t)timeDerivatives) % ALIGNMENT == 0);
  assert(((uintptr_t)timeEvaluated) % ALIGNMENT == 0);

  // assert that this is a forward evaluation in time
  assert(time >= expansionPoint);

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
  }
  intKrnl.power = 1.0;

  // iterate over time derivatives
  for (int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / real(derivative + 1);
  }
}

void seissol::kernels::Time::computeBatchedTaylorExpansion(real time,
                                                           real expansionPoint,
                                                           real** timeDerivatives,
                                                           real** timeEvaluated,
                                                           size_t numElements) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeEvaluated != nullptr);
  assert(time >= expansionPoint);
  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansion::TmpMaxMemRequiredInBytes == 0);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    intKrnl.extraOffset_dQ(i) = m_derivativesOffsets[i];
  }

  // iterate over time derivatives
  const real deltaT = time - expansionPoint;
  intKrnl.power = 1.0;
  for (int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / static_cast<real>(derivative + 1);
  }
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::computeDerivativeTaylorExpansion(real time,
                                                              real expansionPoint,
                                                              real const* timeDerivatives,
                                                              real timeEvaluated[tensor::Q::size()],
                                                              unsigned order) {
  /*
   * assert alignments.
   */
  assert(((uintptr_t)timeDerivatives) % ALIGNMENT == 0);
  assert(((uintptr_t)timeEvaluated) % ALIGNMENT == 0);

  // assert that this is a forward evaluation in time
  assert(time >= expansionPoint);

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
  }
  intKrnl.power = 1.0;

  // iterate over time derivatives
  for (unsigned derivative = order; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / real(derivative + 1);
  }
}

void seissol::kernels::Time::flopsTaylorExpansion(long long& nonZeroFlops,
                                                  long long& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0;
  hardwareFlops = 0;

  // interate over derivatives
  for (unsigned der = 0; der < CONVERGENCE_ORDER; ++der) {
    nonZeroFlops += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
}

unsigned int* seissol::kernels::Time::getDerivativesOffsets() { return m_derivativesOffsets; }

void seissol::kernels::Time::calculateEps(const real *exxNodal, const real *eyyNodal, const real *ezzNodal, const real *exyNodal, const real *eyzNodal, const real *ezxNodal, const unsigned int &q, const seissol::initializer::parameters::DamagedElasticParameters &damagedElasticParameters, real &EspI, real &EspII, real &xi){
      const real epsInitxx = damagedElasticParameters.epsInitxx;
    const real epsInityy = damagedElasticParameters.epsInityy;
    const real epsInitzz = damagedElasticParameters.epsInitzz;
    const real epsInitxy = damagedElasticParameters.epsInitxy;
    const real epsInityz = damagedElasticParameters.epsInitxx;
    const real epsInitzx = damagedElasticParameters.epsInitxx;

    EspI = (exxNodal[q] + epsInitxx) + (eyyNodal[q] + epsInityy) + (ezzNodal[q] + epsInitzz);
    EspII = (exxNodal[q] + epsInitxx) * (exxNodal[q] + epsInitxx) +
                 (eyyNodal[q] + epsInityy) * (eyyNodal[q] + epsInityy) +
                 (ezzNodal[q] + epsInitzz) * (ezzNodal[q] + epsInitzz) +
                 2 * (exyNodal[q] + epsInitxy) * (exyNodal[q] + epsInitxy) +
                 2 * (eyzNodal[q] + epsInityz) * (eyzNodal[q] + epsInityz) +
                 2 * (ezxNodal[q] + epsInitzx) * (ezxNodal[q] + epsInitzx);
    if (EspII > 1e-30) {
      xi = EspI / std::sqrt(EspII);
    } else {
      xi = 0.0;
    }
}

void seissol::kernels::Time::computeNonLinearRusanovFlux(const CellMaterialData *materialData, const unsigned int &l_cell, const unsigned int &side, const double *timeWeights, const real *qIPlus, const real *qIMinus, real *rusanovFluxP, const LocalIntegrationData *localIntegration){
      using namespace seissol::dr::misc::quantity_indices;
    /// Checked that, after reshaping, it still uses the same memory address
    /// S4: Integration in time the Rusanov flux on surface quadrature nodes.
    const unsigned DAM = 9;
    const unsigned BRE = 10;

    const real lambda0P = materialData[l_cell].local.lambda0;
    const real mu0P = materialData[l_cell].local.mu0;
    const real rho0P = materialData[l_cell].local.rho;

    const real lambda0M = materialData[l_cell].neighbor[side].lambda0;
    const real mu0M = materialData[l_cell].neighbor[side].mu0;
    const real rho0M = materialData[l_cell].neighbor[side].rho;

    const real epsInitxx = m_damagedElasticParameters->epsInitxx;
    const real epsInityy = m_damagedElasticParameters->epsInityy;
    const real epsInitzz = m_damagedElasticParameters->epsInitzz;
    const real epsInitxy = m_damagedElasticParameters->epsInitxy;
    const real epsInityz = m_damagedElasticParameters->epsInityz;
    const real epsInitzx = m_damagedElasticParameters->epsInitzx;

    const real aB0 = m_damagedElasticParameters->aB0;
    const real aB1 = m_damagedElasticParameters->aB1;
    const real aB2 = m_damagedElasticParameters->aB2;
    const real aB3 = m_damagedElasticParameters->aB3;

    real lambdaMax = 1.0 * std::sqrt((lambda0P + 2 * mu0P) / rho0P);
    real sxxP, syyP, szzP, sxyP, syzP, szxP, sxxM, syyM, szzM, sxyM, syzM, szxM;

    for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
      auto weight = timeWeights[o];

      for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints; i++) {

        real EspIp, EspIIp, xip, EspIm, EspIIm, xim;

        calculateEps(&qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints],
                                  i,
                                  *m_damagedElasticParameters,
                                  EspIp,
                                  EspIIp,
                                  xip);
        real alphap = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];
        calculateEps(&qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints],
                                  i,
                                  *m_damagedElasticParameters,
                                  EspIm,
                                  EspIIm,
                                  xim);
        real alpham = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];
        real lambp = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                         (lambda0P - alphap * materialData[l_cell].local.gammaR *
                                         (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) / std::sqrt(EspIIp)) +
                     qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (2.0 * aB2 + 3.0 * xip * aB3 +
                                          aB1 * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) / std::sqrt(EspIIp));
        real mup =
            (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                (mu0P -
                 alphap * materialData[l_cell].local.xi0 * materialData[l_cell].local.gammaR -
                 0.5 * alphap * materialData[l_cell].local.gammaR * xip) +
            qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (aB0 + 0.5 * xip * aB1 - 0.5 * xip * xip * xip * aB3);

        real lambm =
            (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                (lambda0M - alpham * materialData[l_cell].neighbor[side].gammaR *
                                (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints+i] + epsInitxx) / std::sqrt(EspIIm)) +
            qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (2.0 * aB2 + 3.0 * xim * aB3 +
                                  aB1 * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) / std::sqrt(EspIIm));

        real mum = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                       (mu0M -
                        alpham * materialData[l_cell].neighbor[side].xi0 *
                            materialData[l_cell].neighbor[side].gammaR -
                        0.5 * alpham * materialData[l_cell].neighbor[side].gammaR * xim) +
                   qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (aB0 + 0.5 * xim * aB1 - 0.5 * xim * xim * xim * aB3);

        lambdaMax =
            std::min(std::sqrt((lambp + 2 * mup) / rho0P), std::sqrt((lambm + 2 * mum) / rho0M));

        // damage stress
        real mu_eff = materialData[l_cell].local.mu0 -
                      alphap * materialData[l_cell].local.gammaR * materialData[l_cell].local.xi0 -
                      0.5 * alphap * materialData[l_cell].local.gammaR * xip;
        real sxx_sp = materialData[l_cell].local.lambda0 * EspIp -
                      alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
                      2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_sp = materialData[l_cell].local.lambda0 * EspIp -
                      alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
                      2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_sp = materialData[l_cell].local.lambda0 * EspIp -
                      alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
                      2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        // breakage stress
        real sxx_bp =
            (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_bp =
            (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_bp =
            (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_bp =
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_bp =
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_bp =
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        // damage stress minus
        mu_eff = materialData[l_cell].neighbor[side].mu0 -
                 alpham * materialData[l_cell].neighbor[side].gammaR *
                     materialData[l_cell].neighbor[side].xi0 -
                 0.5 * alpham * materialData[l_cell].neighbor[side].gammaR * xim;
        real sxx_sm = materialData[l_cell].neighbor[side].lambda0 * EspIm -
                      alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
                      2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_sm = materialData[l_cell].neighbor[side].lambda0 * EspIm -
                      alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
                      2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_sm = materialData[l_cell].neighbor[side].lambda0 * EspIm -
                      alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
                      2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        // breakage stress
        real sxx_bm =
            (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_bm =
            (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_bm =
            (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_bm =
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_bm =
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_bm =
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        real breakp = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i];
        real breakm = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i];

        sxxP = (1 - breakp) * sxx_sp + breakp * sxx_bp;
        syyP = (1 - breakp) * syy_sp + breakp * syy_bp;
        szzP = (1 - breakp) * szz_sp + breakp * szz_bp;
        sxyP = (1 - breakp) * sxy_sp + breakp * sxy_bp;
        syzP = (1 - breakp) * syz_sp + breakp * syz_bp;
        szxP = (1 - breakp) * szx_sp + breakp * szx_bp;

        sxxM = (1 - breakm) * sxx_sm + breakm * sxx_bm;
        syyM = (1 - breakm) * syy_sm + breakm * syy_bm;
        szzM = (1 - breakm) * szz_sm + breakm * szz_bm;

        sxyM = (1 - breakm) * sxy_sm + breakm * sxy_bm;
        syzM = (1 - breakm) * syz_sm + breakm * syz_bm;
        szxM = (1 - breakm) * szx_sm + breakm * szx_bm;

        rusanovFluxP[XX*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx));

        rusanovFluxP[YY*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy));

        rusanovFluxP[ZZ*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz));

        rusanovFluxP[XY*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy));

        rusanovFluxP[YZ*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ( (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz));

        rusanovFluxP[XZ*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx));

        rusanovFluxP[U*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-sxxP / rho0P) + 0.5 * (-sxxM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-sxyP / rho0P) + 0.5 * (-sxyM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-szxP / rho0P) + 0.5 * (-szxM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) - 0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]));

        rusanovFluxP[V*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-sxyP / rho0P) + 0.5 * (-sxyM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-syyP / rho0P) + 0.5 * (-syyM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-syzP / rho0P) + 0.5 * (-syzM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) - 0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]));

        rusanovFluxP[W*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-szxP / rho0P) + 0.5 * (-szxM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-syzP / rho0P) + 0.5 * (-syzM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-szzP / rho0P) + 0.5 * (-szzM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) - 0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]));
      }
    }
}

void seissol::kernels::Time::calculateDynamicRuptureReceiverOutput(const real* dofsNPlus, const seissol::initializer::parameters::DamagedElasticParameters& damagedElasticParameters,
const seissol::dr::ImpedancesAndEta* impAndEtaGet, real* dofsStressNPlus, const real* dofsNMinus, real* dofsStressNMinus){
      for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; q++) {
        real EspIp, EspIIp, EspIm, EspIIm, xip, xim;
        const real epsInitxx = damagedElasticParameters.epsInitxx;
        const real epsInityy = damagedElasticParameters.epsInityy;
        const real epsInitzz = damagedElasticParameters.epsInitzz;
        const real epsInitxy = damagedElasticParameters.epsInitxy;
        const real epsInityz = damagedElasticParameters.epsInityz;
        const real epsInitzx = damagedElasticParameters.epsInitzx;
        const real aB0 = damagedElasticParameters.aB0;
        const real aB1 = damagedElasticParameters.aB1;
        const real aB2 = damagedElasticParameters.aB2;
        const real aB3 = damagedElasticParameters.aB3;

        calculateEps(&dofsNPlus[0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNPlus[1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], &dofsNPlus[2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNPlus[3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], &dofsNPlus[4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNPlus[5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], q, damagedElasticParameters, EspIp, EspIIp, xip);
        real alphap = dofsNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];

        real lambda0P = impAndEtaGet->lambda0P;
        real mu0P = impAndEtaGet->mu0P;
        real lambda0M = impAndEtaGet->lambda0M;
        real mu0M = impAndEtaGet->mu0M;


      // damage stress impAndEtaGet->gammaRP, mu0P
      real mu_eff = mu0P - alphap * impAndEtaGet->gammaRP * impAndEtaGet->xi0P -
                    0.5 * alphap * impAndEtaGet->gammaRP * xip;
      real sxx_sp = lambda0P * EspIp - alphap * impAndEtaGet->gammaRP * std::sqrt(EspIIp) +
                    2 * mu_eff * (dofsNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx);
      real syy_sp = lambda0P * EspIp - alphap * impAndEtaGet->gammaRP * std::sqrt(EspIIp) +
                    2 * mu_eff * (dofsNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy);
      real szz_sp = lambda0P * EspIp - alphap * impAndEtaGet->gammaRP * std::sqrt(EspIIp) +
                    2 * mu_eff * (dofsNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz);

      real sxy_sp = 2 * mu_eff * (dofsNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy);
      real syz_sp = 2 * mu_eff * (dofsNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz);
      real szx_sp = 2 * mu_eff * (dofsNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx);

      // breakage stress
      real sxx_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                    (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
                        (dofsNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx);
      real syy_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                    (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
                        (dofsNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy);
      real szz_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                    (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
                        (dofsNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz);

      real sxy_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
                    (dofsNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy);
      real syz_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
                    (dofsNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz);
      real szx_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
                    (dofsNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx);

      dofsStressNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * sxx_sp +
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * sxx_bp;

      dofsStressNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * syy_sp +
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * syy_bp;

      dofsStressNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * szz_sp +
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * szz_bp;

      dofsStressNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * sxy_sp +
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * sxy_bp;

      dofsStressNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * syz_sp +
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * syz_bp;

      dofsStressNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * szx_sp +
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * szx_bp;

      calculateEps(&dofsNMinus[0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      &dofsNMinus[1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], &dofsNMinus[2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      &dofsNMinus[3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], &dofsNMinus[4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      &dofsNMinus[5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], q, damagedElasticParameters, EspIm, EspIIm, xim);
      real alpham = dofsNMinus[9];

      // damage stress minus
      mu_eff = mu0M - alpham * impAndEtaGet->gammaRM * impAndEtaGet->xi0M -
               0.5 * alpham * impAndEtaGet->gammaRM * xim;
      real sxx_sm =
          lambda0M * EspIm - alpham * impAndEtaGet->gammaRM * std::sqrt(EspIIm) +
          2 * mu_eff * (dofsNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx);
      real syy_sm =
          lambda0M * EspIm - alpham * impAndEtaGet->gammaRM * std::sqrt(EspIIm) +
          2 * mu_eff * (dofsNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy);
      real szz_sm =
          lambda0M * EspIm - alpham * impAndEtaGet->gammaRM * std::sqrt(EspIIm) +
          2 * mu_eff * (dofsNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz);

      real sxy_sm =
          2 * mu_eff * (dofsNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy);
      real syz_sm =
          2 * mu_eff * (dofsNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz);
      real szx_sm =
          2 * mu_eff * (dofsNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx);

      // breakage stress
      real sxx_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                    (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
                        (dofsNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx);
      real syy_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                    (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
                        (dofsNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy);
      real szz_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                    (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
                        (dofsNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz);

      real sxy_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
                    (dofsNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy);
      real syz_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
                    (dofsNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz);
      real szx_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
                    (dofsNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx);

      dofsStressNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * sxx_sm +
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * sxx_bm;

      dofsStressNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * syy_sm +
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * syy_bm;

      dofsStressNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * szz_sm +
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * szz_bm;

      dofsStressNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * sxy_sm +
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * sxy_bm;

      dofsStressNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * syz_sm +
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * syz_bm;

      dofsStressNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (1 - dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q]) * szx_sm +
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] * szx_bm;

      dofsStressNPlus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];

      dofsStressNMinus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
    }
}

void seissol::kernels::Time::computeNonLinearIntegralCorrection(const CellLocalInformation *cellInformation, const unsigned int &l_cell, real **derivatives, real *(*faceNeighbors)[4], const CellMaterialData *materialData, const LocalIntegrationData *localIntegration, const NeighborData &data, const CellDRMapping (*drMapping)[4], kernel::nonlinearSurfaceIntegral &m_nonlSurfIntPrototype, double timeStepSize, const kernel::nonlEvaluateAndRotateQAtInterpolationPoints &m_nonlinearInterpolation){
  double timePoints[CONVERGENCE_ORDER];
      double timeWeights[CONVERGENCE_ORDER];

      seissol::quadrature::GaussLegendre(timePoints, timeWeights, CONVERGENCE_ORDER);

      for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
        timePoints[point] = 0.5 * (timeStepSize * timePoints[point] + timeStepSize);
        timeWeights[point] = 0.5 * timeStepSize * timeWeights[point];
      }
        for (unsigned int side = 0; side < 4; side++ ){
          if (cellInformation[l_cell].faceTypes[side] == FaceType::regular
          || cellInformation[l_cell].faceTypes[side] == FaceType::periodic){
            // Compute local integrals with derivatives and Rusanov flux
            /// S1: compute the space-time interpolated Q on both side of 4 faces
            /// S2: at the same time rotate the field to face-aligned coord.
            alignas(PAGESIZE_STACK) real QInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};
            alignas(PAGESIZE_STACK) real QInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};

            for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
              real degreesOfFreedomPlus[tensor::Q::size()];
              real degreesOfFreedomMinus[tensor::Q::size()];

              for (unsigned i_f = 0; i_f < tensor::Q::size(); i_f++){
                degreesOfFreedomPlus[i_f] = static_cast<real>(0.0);
                degreesOfFreedomMinus[i_f] = static_cast<real>(0.0);
              }

              // !!! Make sure every time after entering this function, the last input should be reinitialized to zero
              computeTaylorExpansion(timePoints[timeInterval], 0.0, derivatives[l_cell], degreesOfFreedomPlus);
              computeTaylorExpansion(timePoints[timeInterval], 0.0, faceNeighbors[l_cell][side], degreesOfFreedomMinus);

              // Prototype is necessary for openmp
              kernel::nonlEvaluateAndRotateQAtInterpolationPoints m_nonLinInter
                = m_nonlinearInterpolation;

              m_nonLinInter.QInterpolated = &QInterpolatedPlus[timeInterval][0];
              m_nonLinInter.Q = degreesOfFreedomPlus;
              m_nonLinInter.execute(side, 0);

              m_nonLinInter.QInterpolated = &QInterpolatedMinus[timeInterval][0];
              m_nonLinInter.Q = degreesOfFreedomMinus;
              m_nonLinInter.execute(cellInformation[l_cell].faceRelations[side][0]
                            , cellInformation[l_cell].faceRelations[side][1]+1);
            }

            // S3: Construct matrices to store Rusanov flux on surface quadrature nodes.
            // Reshape the interpolated results
            using QInterpolatedShapeT = const real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];

            auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedPlus));
            auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedMinus));

            // The arrays to store time integrated flux
            alignas(PAGESIZE_STACK) real rusanovFluxPlus[tensor::QInterpolated::size()] = {0.0};
            // alignas(PAGESIZE_STACK) real rusanovFluxMinus[tensor::QInterpolated::size()] = {0.0};

            for (unsigned i_f = 0; i_f < tensor::QInterpolated::size(); i_f++){
              rusanovFluxPlus[i_f] = static_cast<real>(0.0);
            }

            using rusanovFluxShape = real(*)[seissol::dr::misc::numPaddedPoints];
            auto* rusanovFluxP = reinterpret_cast<rusanovFluxShape>(rusanovFluxPlus);

            // S4: Compute the Rusanov flux
            computeNonLinearRusanovFlux(materialData, l_cell, side, timeWeights, *qIPlus[0], *qIMinus[0], *rusanovFluxP, localIntegration);

            /// S5: Integrate in space using quadrature.
            kernel::nonlinearSurfaceIntegral m_surfIntegral = m_nonlSurfIntPrototype;
            m_surfIntegral.Q = data.dofs;
            m_surfIntegral.Flux = rusanovFluxPlus;
            m_surfIntegral.fluxScale = localIntegration[l_cell].fluxScales[side];
            m_surfIntegral.execute(side, 0);
          }
          else if (cellInformation[l_cell].faceTypes[side] == FaceType::dynamicRupture) {
            // No neighboring cell contribution, interior bc.
            assert(reinterpret_cast<uintptr_t>(drMapping[l_cell][side].godunov) % ALIGNMENT == 0);

            kernel::nonlinearSurfaceIntegral m_drIntegral = m_nonlSurfIntPrototype;
            m_drIntegral.Q = data.dofs;
            m_drIntegral.Flux = drMapping[l_cell][side].godunov;
            m_drIntegral.fluxScale = localIntegration[l_cell].fluxScales[side];
            m_drIntegral.execute(side, drMapping[l_cell][side].faceRelation);
          } // if (faceTypes)
        } // for (side)
}

void seissol::kernels::Time::stressToDofsDynamicRupture(real* dofsStressNPlus, const real* dofsNPlus, real* dofsStressNMinus, 
const real* dofsNMinus){
    const real epsInitxx = m_damagedElasticParameters->epsInitxx;
    const real epsInityy = m_damagedElasticParameters->epsInityy;
    const real epsInitzz = m_damagedElasticParameters->epsInitzz;
    const real epsInitxy = m_damagedElasticParameters->epsInitxy;
    const real epsInityz = m_damagedElasticParameters->epsInityz;
    const real epsInitzx = m_damagedElasticParameters->epsInitzx;

    for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; q++) {
      dofsStressNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx);

      dofsStressNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy);

      dofsStressNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz);

      dofsStressNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy);

      dofsStressNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz);

      dofsStressNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx);

      dofsStressNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx);

      dofsStressNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy);

      dofsStressNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz);

      dofsStressNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy);

      dofsStressNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz);

      dofsStressNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          (dofsNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx);

      dofsStressNPlus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];

      dofsStressNMinus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
      dofsStressNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
    }
}

void seissol::kernels::Time::computeNonLinearBaseFrictionLaw(const seissol::dr::ImpedancesAndEta* impAndEta, const unsigned& ltsFace,
const real* qIPlus, real* qStressIPlus, const real* qIMinus, real* qStressIMinus){
      using namespace seissol::dr::misc::quantity_indices;
      const real lambda0P = impAndEta[ltsFace].lambda0P;
      const real mu0P = impAndEta[ltsFace].mu0P;
      const real lambda0M = impAndEta[ltsFace].lambda0M;
      const real mu0M = impAndEta[ltsFace].mu0M;

      // TODO(NONLINEAR) What are these values?
      const real aB0 = m_damagedElasticParameters->aB0;
      const real aB1 = m_damagedElasticParameters->aB1;
      const real aB2 = m_damagedElasticParameters->aB2;
      const real aB3 = m_damagedElasticParameters->aB3;

      for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
        for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints; i++) {
          real EspIp, EspIIp, alphap, xip, EspIm, EspIIm, alpham, xim;
          calculateEps(&qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints],
                                  i,
                                  *m_damagedElasticParameters,
                                  EspIp,
                                  EspIIp,
                                  xip);
          alphap = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];
          // damage stress impAndEtaGet->gammaRP, mu0P
          real mu_eff = mu0P - alphap * impAndEta[ltsFace].gammaRP * impAndEta[ltsFace].xi0P -
                        0.5 * alphap * impAndEta[ltsFace].gammaRP * xip;
          real sxx_sp = lambda0P * EspIp - alphap * impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp) +
                        2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i]);
          real syy_sp = lambda0P * EspIp - alphap * impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp) +
                        2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i]);
          real szz_sp = lambda0P * EspIp - alphap * impAndEta[ltsFace].gammaRP * std::sqrt(EspIIp) +
                        2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i]);

          const real sxy_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i]);
          const real syz_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i]);
          const real szx_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i]);

          // breakage stress
          const real sxx_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                              (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i]);
          const real syy_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                              (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i]);
          const real szz_bp = (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
                              (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i]);

          const real sxy_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i]);
          const real syz_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i]);
          const real szx_bp = (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i]);

          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * sxx_sp + qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * sxx_bp;
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * syy_sp + qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * syy_bp;
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * szz_sp + qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * szz_bp;
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * sxy_sp + qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * sxy_bp;
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * syz_sp + qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * syz_bp;
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * szx_sp + qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * szx_bp;

          calculateEps(&qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints],
                                  i,
                                  *m_damagedElasticParameters,
                                  EspIm,
                                  EspIIm,
                                  xim);
          alpham = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];

          // damage stress minus
          mu_eff = mu0M - alpham * impAndEta[ltsFace].gammaRM * impAndEta[ltsFace].xi0M -
                   0.5 * alpham * impAndEta[ltsFace].gammaRM * xim;
          const real sxx_sm = lambda0M * EspIm -
                              alpham * impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm) +
                              2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i]);
          const real syy_sm = lambda0M * EspIm -
                              alpham * impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm) +
                              2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i]);
          const real szz_sm = lambda0M * EspIm -
                              alpham * impAndEta[ltsFace].gammaRM * std::sqrt(EspIIm) +
                              2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i]);

          const real sxy_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i]);
          const real syz_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i]);
          const real szx_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i]);

          // breakage stress
          const real sxx_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                              (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i]);
          const real syy_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                              (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i]);
          const real szz_bm = (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
                              (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i]);

          const real sxy_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i]);
          const real syz_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i]);
          const real szx_bm = (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i]);


          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * sxx_sm + qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * sxx_bm;
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * syy_sm + qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * syy_bm;
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * szz_sm + qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * szz_bm;
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * sxy_sm + qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * sxy_bm;
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * syz_sm + qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * syz_bm;
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) * szx_sm + qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * szx_bm;

          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i] = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i];
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i] = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i];
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i] = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i];
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i] = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];
          qStressIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i];

          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i] = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i];
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i] = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i];
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i] = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i];
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i] = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];
          qStressIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i];
        }
      }
}