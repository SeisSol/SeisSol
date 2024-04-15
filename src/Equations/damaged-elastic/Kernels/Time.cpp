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
#include <Numerical_aux/Transformation.h>
#include <tuple>

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
    std::tie(EspI, EspII, xi) = calculateEsp(exxNodal,
                 eyyNodal,
                 ezzNodal,
                 exyNodal,
                 eyzNodal,
                 ezxNodal,
                 q,
                 m_damagedElasticParameters);

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

std::tuple<real, real, real> seissol::kernels::Time::calculateEsp(
    const real* exxNodal,
    const real* eyyNodal,
    const real* ezzNodal,
    const real* exyNodal,
    const real* eyzNodal,
    const real* ezxNodal,
    unsigned int q,
    const seissol::initializer::parameters::DamagedElasticParameters* damagedElasticParameters) {
  const real epsInitxx = damagedElasticParameters->epsInitxx;
  const real epsInityy = damagedElasticParameters->epsInityy;
  const real epsInitzz = damagedElasticParameters->epsInitzz;
  const real epsInitxy = damagedElasticParameters->epsInitxy;
  const real epsInityz = damagedElasticParameters->epsInitxx;
  const real epsInitzx = damagedElasticParameters->epsInitxx;

  real EspI = (exxNodal[q] + epsInitxx) + (eyyNodal[q] + epsInityy) + (ezzNodal[q] + epsInitzz);
  real EspII = (exxNodal[q] + epsInitxx) * (exxNodal[q] + epsInitxx) +
          (eyyNodal[q] + epsInityy) * (eyyNodal[q] + epsInityy) +
          (ezzNodal[q] + epsInitzz) * (ezzNodal[q] + epsInitzz) +
          2 * (exyNodal[q] + epsInitxy) * (exyNodal[q] + epsInitxy) +
          2 * (eyzNodal[q] + epsInityz) * (eyzNodal[q] + epsInityz) +
          2 * (ezxNodal[q] + epsInitzx) * (ezxNodal[q] + epsInitzx);
    real xi;
  if (EspII > 1e-30) {
    xi = EspI / std::sqrt(EspII);
  } else {
    xi = 0.0;
  }
  return {EspI, EspII, xi};
}

void seissol::kernels::Time::calculateDynamicRuptureReceiverOutput(
    const real* dofsNPlus,
    const seissol::initializer::parameters::DamagedElasticParameters& damagedElasticParameters,
    const seissol::dr::ImpedancesAndEta* impAndEtaGet,
    real* dofsStressNPlus,
    const real* dofsNMinus,
    real* dofsStressNMinus) {
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

    std::tie(EspIp, EspIIp, xip) = calculateEsp(&dofsNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 q,
                 &damagedElasticParameters);
    real alphap = dofsNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];

    real lambda0P = impAndEtaGet->lambda0P;
    real mu0P = impAndEtaGet->mu0P;
    real lambda0M = impAndEtaGet->lambda0M;
    real mu0M = impAndEtaGet->mu0M;

    real mu_eff, sxx_sp, syy_sp, szz_sp, sxy_sp, syz_sp, szx_sp, sxx_bp, syy_bp, szz_bp, sxy_bp,
        syz_bp, szx_bp;

    std::tie(mu_eff,
             sxx_sp,
             syy_sp,
             szz_sp,
             sxy_sp,
             syz_sp,
             szx_sp,
             sxx_bp,
             syy_bp,
             szz_bp,
             sxy_bp,
             syz_bp,
             szx_bp) =
        calculateDamageAndBreakageStresses(
            mu0P,
            alphap,
            impAndEtaGet->gammaRP,
            impAndEtaGet->xi0P,
            xip,
            lambda0P,
            EspIp,
            EspIIp,
            dofsNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx,
            dofsNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy,
            dofsNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz,
            dofsNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy,
            dofsNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz,
            dofsNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx,
            aB0,
            aB1,
            aB2,
            aB3);

calculateStressesFromDamageAndBreakageStresses(&dofsStressNPlus[0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    &dofsStressNPlus[3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    &dofsStressNPlus[4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[6*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[7*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[8*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNPlus[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNPlus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    &dofsNPlus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNPlus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    sxx_sp, sxx_bp, syy_sp, syy_bp, szz_sp, szz_bp, sxy_sp, sxy_bp, syz_sp, syz_bp, szx_sp, szx_bp, q);

    std::tie(EspIm, EspIIm, xim) = calculateEsp(&dofsNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 &dofsNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                 q,
                 &damagedElasticParameters);
    real alpham = dofsNMinus[9];

    real sxx_sm, syy_sm, szz_sm, sxy_sm, syz_sm, szx_sm, sxx_bm, syy_bm, szz_bm, sxy_bm, syz_bm,
        szx_bm;
    std::tie(mu_eff,
             sxx_sm,
             syy_sm,
             szz_sm,
             sxy_sm,
             syz_sm,
             szx_sm,
             sxx_bm,
             syy_bm,
             szz_bm,
             sxy_bm,
             syz_bm,
             szx_bm) =
        calculateDamageAndBreakageStresses(
            mu0M,
            alpham,
            impAndEtaGet->gammaRM,
            impAndEtaGet->xi0M,
            xim,
            lambda0M,
            EspIm,
            EspIIm,
            dofsNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx,
            dofsNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy,
            dofsNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz,
            dofsNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy,
            dofsNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz,
            dofsNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx,
            aB0,
            aB1,
            aB2,
            aB3);

calculateStressesFromDamageAndBreakageStresses(&dofsStressNMinus[0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    &dofsStressNMinus[3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    &dofsStressNMinus[4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[6*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[7*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[8*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsStressNMinus[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNMinus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    &dofsNMinus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNMinus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNMinus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS], 
    &dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    sxx_sm, sxx_bm, syy_sm, syy_bm, szz_sm, szz_bm, sxy_sm, sxy_bm, syz_sm, syz_bm, szx_sm, szx_bm, q);
  }
}

void seissol::kernels::Time::stressToDofsDynamicRupture(real* dofsStress,
                                                        const real* dofs) {
  const real epsInitxx = m_damagedElasticParameters->epsInitxx;
  const real epsInityy = m_damagedElasticParameters->epsInityy;
  const real epsInitzz = m_damagedElasticParameters->epsInitzz;
  const real epsInitxy = m_damagedElasticParameters->epsInitxy;
  const real epsInityz = m_damagedElasticParameters->epsInityz;
  const real epsInitzx = m_damagedElasticParameters->epsInitzx;

  for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; q++) {
    dofsStress[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        (dofs[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxx);

    dofsStress[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        (dofs[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityy);

    dofsStress[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        (dofs[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzz);

    dofsStress[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        (dofs[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitxy);

    dofsStress[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        (dofs[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInityz);

    dofsStress[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        (dofs[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] + epsInitzx);

    dofsStress[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        dofs[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
    dofsStress[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        dofs[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
    dofsStress[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        dofs[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
    dofsStress[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        dofs[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
    dofsStress[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        dofs[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q];
  }
}

void seissol::kernels::Time::computeNonLinearBaseFrictionLaw(
    const seissol::dr::ImpedancesAndEta* impAndEta,
    unsigned ltsFace,
    const real* qIPlus,
    real* qStressIPlus,
    const real* qIMinus,
    real* qStressIMinus) {
  using namespace seissol::dr::misc::quantity_indices;
  using namespace seissol::dr::misc;

  real lambda0P = impAndEta[ltsFace].lambda0P;
  real mu0P = impAndEta[ltsFace].mu0P;
  real lambda0M = impAndEta[ltsFace].lambda0M;
  real mu0M = impAndEta[ltsFace].mu0M;

  // TODO(NONLINEAR) What are these values?
  const real aB0 = m_damagedElasticParameters->aB0;
  const real aB1 = m_damagedElasticParameters->aB1;
  const real aB2 = m_damagedElasticParameters->aB2;
  const real aB3 = m_damagedElasticParameters->aB3;

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints; i++) {

      real EspIp = (qIPlus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i]) +
                   (qIPlus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i]) +
                   (qIPlus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i]);
      real EspIIp = (qIPlus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i]) *
                        (qIPlus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i]) +
                    (qIPlus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i]) *
                        (qIPlus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i]) +
                    (qIPlus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i]) *
                        (qIPlus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i]) +
                    2 * (qIPlus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints + i]) *
                        (qIPlus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints + i]) +
                    2 * (qIPlus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints + i]) *
                        (qIPlus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints + i]) +
                    2 * (qIPlus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints + i]) *
                        (qIPlus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints + i]);
      real alphap = qIPlus[o * numQuantities * numPaddedPoints + DAM * numPaddedPoints + i];
      real xip;
      if (EspIIp > 1e-30) {
        xip = EspIp / std::sqrt(EspIIp);
      } else {
        xip = 0.0;
      }

    real mu_eff, sxx_sp, syy_sp, szz_sp, sxy_sp, syz_sp, szx_sp, sxx_bp, syy_bp, szz_bp, sxy_bp,
        syz_bp, szx_bp;
    std::tie(mu_eff,
             sxx_sp,
             syy_sp,
             szz_sp,
             sxy_sp,
             syz_sp,
             szx_sp,
             sxx_bp,
             syy_bp,
             szz_bp,
             sxy_bp,
             syz_bp,
             szx_bp) =
        calculateDamageAndBreakageStresses(
            mu0P,
            alphap,
            impAndEta[ltsFace].gammaRP,
            impAndEta[ltsFace].xi0P,
            xip,
            lambda0P,
            EspIp,
            EspIIp,
            qIPlus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i],
            qIPlus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i],
            qIPlus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i],
            qIPlus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints + i],
            qIPlus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints + i],
            qIPlus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints + i],
            aB0,
            aB1,
            aB2,
            aB3);
    
    calculateStressesFromDamageAndBreakageStresses(&qStressIPlus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + U * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + V * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + W * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + DAM * numPaddedPoints], 
    &qStressIPlus[o * numQuantities * numPaddedPoints + BRE * numPaddedPoints], 
    &qIPlus[o * numQuantities * numPaddedPoints + U * numPaddedPoints], 
    &qIPlus[o * numQuantities * numPaddedPoints + V * numPaddedPoints], 
    &qIPlus[o * numQuantities * numPaddedPoints + W * numPaddedPoints], 
    &qIPlus[o * numQuantities * numPaddedPoints + DAM * numPaddedPoints], 
    &qIPlus[o * numQuantities * numPaddedPoints + BRE * numPaddedPoints], 
    sxx_sp, sxx_bp, syy_sp, syy_bp, szz_sp, szz_bp, sxy_sp, sxy_bp, syz_sp, syz_bp, szx_sp, szx_bp, i);

      const real EspIm = (qIMinus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i]) +
                         (qIMinus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i]) +
                         (qIMinus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i]);
      const real EspIIm =
          (qIMinus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i]) *
              (qIMinus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i]) +
          (qIMinus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i]) *
              (qIMinus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i]) +
          (qIMinus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i]) *
              (qIMinus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i]) +
          2 * (qIMinus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints + i]) *
              (qIMinus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints + i]) +
          2 * (qIMinus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints + i]) *
              (qIMinus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints + i]) +
          2 * (qIMinus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints + i]) *
              (qIMinus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints + i]);
      real alpham = qIMinus[o * numQuantities * numPaddedPoints + DAM * numPaddedPoints + i];
      real xim;
      if (EspIIm > 1e-30) {
        xim = EspIm / std::sqrt(EspIIm);
      } else {
        xim = 0.0;
      }

    real sxx_sm, syy_sm, szz_sm, sxy_sm, syz_sm, szx_sm, sxx_bm, syy_bm, szz_bm, sxy_bm,
        syz_bm, szx_bm;

    std::tie(mu_eff,
             sxx_sm,
             syy_sm,
             szz_sm,
             sxy_sm,
             syz_sm,
             szx_sm,
             sxx_bm,
             syy_bm,
             szz_bm,
             sxy_bm,
             syz_bm,
             szx_bm) =
        calculateDamageAndBreakageStresses(
            mu0M,
            alpham,
            impAndEta[ltsFace].gammaRM,
            impAndEta[ltsFace].xi0M,
            xim,
            lambda0M,
            EspIm,
            EspIIm,
            qIMinus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints + i],
            qIMinus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints + i],
            qIMinus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints + i],
            qIMinus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints + i],
            qIMinus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints + i],
            qIMinus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints + i],
            aB0,
            aB1,
            aB2,
            aB3);

    calculateStressesFromDamageAndBreakageStresses(&qStressIMinus[o * numQuantities * numPaddedPoints + XX * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + YY * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + ZZ * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + XY * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + YZ * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + XZ * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + U * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + V * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + W * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + DAM * numPaddedPoints], 
    &qStressIMinus[o * numQuantities * numPaddedPoints + BRE * numPaddedPoints], 
    &qIMinus[o * numQuantities * numPaddedPoints + U * numPaddedPoints], 
    &qIMinus[o * numQuantities * numPaddedPoints + V * numPaddedPoints], 
    &qIMinus[o * numQuantities * numPaddedPoints + W * numPaddedPoints], 
    &qIMinus[o * numQuantities * numPaddedPoints + DAM * numPaddedPoints], 
    &qIMinus[o * numQuantities * numPaddedPoints + BRE * numPaddedPoints], 
    sxx_sm, sxx_bm, syy_sm, syy_bm, szz_sm, szz_bm, sxy_sm, sxy_bm, syz_sm, syz_bm, szx_sm, szx_bm, i);
    }
  } // time integration loop
}

void seissol::kernels::Time::updateNonLinearMaterialLocal(const real* Q_aveData,
                                                          CellMaterialData* materialData,
                                                          unsigned int l_cell,
                                                          const std::vector<Element>& elements,
                                                          const std::vector<Vertex>& vertices,
                                                          unsigned int meshId,
                                                          real* x,
                                                          real* y,
                                                          real* z) {
  const real epsInitxx = m_damagedElasticParameters->epsInitxx;
  const real epsInityy = m_damagedElasticParameters->epsInityy;
  const real epsInitzz = m_damagedElasticParameters->epsInitzz;
  const real epsInitxy = m_damagedElasticParameters->epsInitxy;
  const real epsInityz = m_damagedElasticParameters->epsInityz;
  const real epsInitzx = m_damagedElasticParameters->epsInitzx;

  const real EspI =
      (Q_aveData[0] + epsInitxx) + (Q_aveData[1] + epsInityy) + (Q_aveData[2] + epsInitzz);
  const real EspII = (Q_aveData[0] + epsInitxx) * (Q_aveData[0] + epsInitxx) +
                     (Q_aveData[1] + epsInityy) * (Q_aveData[1] + epsInityy) +
                     (Q_aveData[2] + epsInitzz) * (Q_aveData[2] + epsInitzz) +
                     2 * (Q_aveData[3] + epsInitxy) * (Q_aveData[3] + epsInitxy) +
                     2 * (Q_aveData[4] + epsInityz) * (Q_aveData[4] + epsInityz) +
                     2 * (Q_aveData[5] + epsInitzx) * (Q_aveData[5] + epsInitzx);

  real xi;
  if (EspII > 1e-30) {
    xi = EspI / std::sqrt(EspII);
  } else {
    xi = 0.0;
  }

  const real alphaAve = Q_aveData[9];
  const real breakAve = Q_aveData[10];

  real lambda0 = materialData[l_cell].local.lambda0;
  real mu0 = materialData[l_cell].local.mu0;

  const real aB0 = m_damagedElasticParameters->aB0;
  const real aB1 = m_damagedElasticParameters->aB1;
  const real aB2 = m_damagedElasticParameters->aB2;
  const real aB3 = m_damagedElasticParameters->aB3;

  materialData[l_cell].local.mu =
      (1 - breakAve) *
          (mu0 - alphaAve * materialData[l_cell].local.xi0 * materialData[l_cell].local.gammaR -
           0.5 * alphaAve * materialData[l_cell].local.gammaR * xi) +
      breakAve * ((aB0 + 0.5 * aB1 * xi - 0.5 * aB3 * xi * xi * xi));
  materialData[l_cell].local.lambda =
      (1 - breakAve) * (lambda0 - alphaAve * materialData[l_cell].local.gammaR *
                                      (Q_aveData[0] + epsInitxx) / std::sqrt(EspII)) +
      breakAve *
          ((2.0 * aB2 + 3.0 * aB3 * xi) + aB1 * (Q_aveData[0] + epsInitxx) / std::sqrt(EspII));
  materialData[l_cell].local.gamma = alphaAve * materialData[l_cell].local.gammaR;

  materialData[l_cell].local.epsxx_alpha = (Q_aveData[0] + epsInitxx);
  materialData[l_cell].local.epsyy_alpha = (Q_aveData[1] + epsInityy);
  materialData[l_cell].local.epszz_alpha = (Q_aveData[2] + epsInitzz);
  materialData[l_cell].local.epsxy_alpha = (Q_aveData[3] + epsInitxy);
  materialData[l_cell].local.epsyz_alpha = (Q_aveData[4] + epsInityz);
  materialData[l_cell].local.epszx_alpha = (Q_aveData[5] + epsInitzx);
}

void seissol::kernels::Time::updateNonLinearMaterialNeighbor(CellMaterialData* materialData,
                                                             unsigned int l_cell,
                                                             unsigned side,
                                                             const real* Q_aveData) {
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

  real lambda0 = materialData[l_cell].neighbor[side].lambda0;
  real mu0 = materialData[l_cell].neighbor[side].mu0;

  real EspINeigh =
      ((Q_aveData[0] + epsInitxx) + (Q_aveData[1] + epsInityy) + (Q_aveData[2] + epsInitzz));
  real EspII = (Q_aveData[0] + epsInitxx) * (Q_aveData[0] + epsInitxx) +
               (Q_aveData[1] + epsInityy) * (Q_aveData[1] + epsInityy) +
               (Q_aveData[2] + epsInitzz) * (Q_aveData[2] + epsInitzz) +
               2 * (Q_aveData[3] + epsInitxy) * (Q_aveData[3] + epsInitxy) +
               2 * (Q_aveData[4] + epsInityz) * (Q_aveData[4] + epsInityz) +
               2 * (Q_aveData[5] + epsInitzx) * (Q_aveData[5] + epsInitzx);
  const real alphaAveNeigh = Q_aveData[9];
  const real breakAveNeigh = Q_aveData[10];
  real xi;
  if (EspII > 1e-30) {
    xi = EspINeigh / std::sqrt(EspII);
  } else {
    xi = 0.0;
  }

  materialData[l_cell].neighbor[side].mu =
      (1 - breakAveNeigh) *
          (mu0 -
           alphaAveNeigh * materialData[l_cell].neighbor[side].xi0 *
               materialData[l_cell].neighbor[side].gammaR -
           0.5 * alphaAveNeigh * materialData[l_cell].neighbor[side].gammaR * xi) +
      breakAveNeigh * ((aB0 + 0.5 * aB1 * xi - 0.5 * aB3 * xi * xi * xi));
  materialData[l_cell].neighbor[side].lambda =
      (1 - breakAveNeigh) * (lambda0 - alphaAveNeigh * materialData[l_cell].neighbor[side].gammaR *
                                           (Q_aveData[0] + epsInitxx) / std::sqrt(EspII)) +
      breakAveNeigh *
          ((2.0 * aB2 + 3.0 * aB3 * xi) + aB1 * (Q_aveData[0] + epsInitxx) / std::sqrt(EspII));
  materialData[l_cell].neighbor[side].gamma =
      alphaAveNeigh * materialData[l_cell].neighbor[side].gammaR;

  materialData[l_cell].neighbor[side].epsxx_alpha = (Q_aveData[0] + epsInitxx);
  materialData[l_cell].neighbor[side].epsyy_alpha = (Q_aveData[1] + epsInityy);
  materialData[l_cell].neighbor[side].epszz_alpha = (Q_aveData[2] + epsInitzz);
  materialData[l_cell].neighbor[side].epsxy_alpha = (Q_aveData[3] + epsInitxy);
  materialData[l_cell].neighbor[side].epsyz_alpha = (Q_aveData[4] + epsInityz);
  materialData[l_cell].neighbor[side].epszx_alpha = (Q_aveData[5] + epsInitzx);
}

void seissol::kernels::Time::computeNonLinearLocalIntegration(
    const seissol::kernels::LocalData& data,
    real (&QInterpolatedBodyNodal)[CONVERGENCE_ORDER][tensor::QNodal::size()],
    real (&FInterpolatedBody)[CONVERGENCE_ORDER][tensor::QNodal::size()],
    real (&sxxNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    real (&syyNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    real (&szzNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    real (&sxyNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    real (&syzNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    real (&szxNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
    real (&FluxInterpolatedBodyX)[CONVERGENCE_ORDER][tensor::QNodal::size()],
    real (&FluxInterpolatedBodyY)[CONVERGENCE_ORDER][tensor::QNodal::size()],
    real (&FluxInterpolatedBodyZ)[CONVERGENCE_ORDER][tensor::QNodal::size()]) {
  const real epsInitxx = m_damagedElasticParameters->epsInitxx;
  const real epsInityy = m_damagedElasticParameters->epsInityy;
  const real epsInitzz = m_damagedElasticParameters->epsInitzz;
  const real epsInitxy = m_damagedElasticParameters->epsInitxy;
  const real epsInityz = m_damagedElasticParameters->epsInityz;
  const real epsInitzx = m_damagedElasticParameters->epsInitzx;
  const real damageParameter = data.material.local.Cd;
  const real scalingValue = m_damagedElasticParameters->scalingValue;
  const real breakCoefficient = scalingValue * damageParameter;
  const real betaAlpha = m_damagedElasticParameters->betaAlpha;
  const real aB0 = m_damagedElasticParameters->aB0;
  const real aB1 = m_damagedElasticParameters->aB1;
  const real aB2 = m_damagedElasticParameters->aB2;
  const real aB3 = m_damagedElasticParameters->aB3;

  for (unsigned int timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    real* exxNodal = (QInterpolatedBodyNodal[timeInterval] + 0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* eyyNodal = (QInterpolatedBodyNodal[timeInterval] + 1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* ezzNodal = (QInterpolatedBodyNodal[timeInterval] + 2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* alphaNodal =
        (QInterpolatedBodyNodal[timeInterval] + 9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* breakNodal =
        (QInterpolatedBodyNodal[timeInterval] + 10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);

    real* exyNodal = (QInterpolatedBodyNodal[timeInterval] + 3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* eyzNodal = (QInterpolatedBodyNodal[timeInterval] + 4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* ezxNodal = (QInterpolatedBodyNodal[timeInterval] + 5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* vxNodal = (QInterpolatedBodyNodal[timeInterval] + 6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* vyNodal = (QInterpolatedBodyNodal[timeInterval] + 7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
    real* vzNodal = (QInterpolatedBodyNodal[timeInterval] + 8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);

    real alpha_ave = alphaNodal[0];
    real break_ave = breakNodal[0];
    for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS - 1; ++q) {
      break_ave = std::max(break_ave, breakNodal[q]);
      alpha_ave = std::max(alpha_ave, alphaNodal[q]);
    }

    for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q) {
      real EspI = (exxNodal[q] + epsInitxx) + (eyyNodal[q] + epsInityy) + (ezzNodal[q] + epsInitzz);
      real EspII = (exxNodal[q] + epsInitxx) * (exxNodal[q] + epsInitxx) +
                   (eyyNodal[q] + epsInityy) * (eyyNodal[q] + epsInityy) +
                   (ezzNodal[q] + epsInitzz) * (ezzNodal[q] + epsInitzz) +
                   2 * (exyNodal[q] + epsInitxy) * (exyNodal[q] + epsInitxy) +
                   2 * (eyzNodal[q] + epsInityz) * (eyzNodal[q] + epsInityz) +
                   2 * (ezxNodal[q] + epsInitzx) * (ezxNodal[q] + epsInitzx);
      real xi;
      if (EspII > 1e-30) {
        xi = EspI / std::sqrt(EspII);
      } else {
        xi = 0.0;
      }

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

      const real alphaCR1q = (-bCR - std::sqrt(bCR * bCR - 4.0 * aCR * cCR)) / (2.0 * aCR);
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
            FInterpolatedBody[timeInterval][10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
                (1 - breakNodal[q]) * 1.0 /
                (std::exp((alphaCRq - alphaNodal[q]) / betaAlpha) + 1.0) * breakCoefficient *
                data.material.local.gammaR * EspII * (xi + data.material.local.xi0);
            FInterpolatedBody[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
                (1 - breakNodal[q]) * damageParameter * data.material.local.gammaR * EspII *
                (xi + data.material.local.xi0);
          } else {
            FInterpolatedBody[timeInterval][10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
            FInterpolatedBody[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
          }
        } else {
          FInterpolatedBody[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
          FInterpolatedBody[timeInterval][10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
        }
      } else if (alpha_ave > 5e-1) {
        FInterpolatedBody[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
            0.0 * damageParameter * data.material.local.gammaR * EspII *
            (xi + data.material.local.xi0);
        FInterpolatedBody[timeInterval][10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
            0.0 * damageParameter * data.material.local.gammaR * EspII *
            (xi + data.material.local.xi0);
      } else {
        FInterpolatedBody[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
        FInterpolatedBody[timeInterval][10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      }

    real mu_eff, sxx_s, syy_s, szz_s, sxy_s, syz_s, szx_s, sxx_b, syy_b, szz_b, sxy_b,
        syz_b, szx_b;
    std::tie(mu_eff,
             sxx_s,
             syy_s,
             szz_s,
             sxy_s,
             syz_s,
             szx_s,
             sxx_b,
             syy_b,
             szz_b,
             sxy_b,
             syz_b,
             szx_b) =
        calculateDamageAndBreakageStresses(
            data.material.local.mu0,
            alphaNodal[q],
            data.material.local.gammaR,
            data.material.local.xi0,
            xi,
            data.material.local.lambda0,
            EspI,
            EspII,
            exxNodal[q]+epsInitxx,
            eyyNodal[q]+epsInityy,
            ezzNodal[q]+epsInitzz,
            exyNodal[q]+epsInitxy,
            eyzNodal[q]+epsInityz,
            ezxNodal[q]+epsInitzx,
            aB0,
            aB1,
            aB2,
            aB3);

      sxxNodal[q] = (1 - breakNodal[q]) * sxx_s + breakNodal[q] * sxx_b;
      syyNodal[q] = (1 - breakNodal[q]) * syy_s + breakNodal[q] * syy_b;
      szzNodal[q] = (1 - breakNodal[q]) * szz_s + breakNodal[q] * szz_b;
      sxyNodal[q] = (1 - breakNodal[q]) * sxy_s + breakNodal[q] * sxy_b;
      syzNodal[q] = (1 - breakNodal[q]) * syz_s + breakNodal[q] * syz_b;
      szxNodal[q] = (1 - breakNodal[q]) * szx_s + breakNodal[q] * szx_b;

      // //--- x-dir
      FluxInterpolatedBodyX[timeInterval][0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -vxNodal[q];
      FluxInterpolatedBodyX[timeInterval][1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyX[timeInterval][2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyX[timeInterval][3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -0.5 * vyNodal[q];
      FluxInterpolatedBodyX[timeInterval][4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyX[timeInterval][5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -0.5 * vzNodal[q];
      FluxInterpolatedBodyX[timeInterval][6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -sxxNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyX[timeInterval][7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -sxyNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyX[timeInterval][8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -szxNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyX[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      //--- y-dir
      FluxInterpolatedBodyY[timeInterval][0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyY[timeInterval][1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -vyNodal[q];
      FluxInterpolatedBodyY[timeInterval][2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyY[timeInterval][3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -0.5 * vxNodal[q];
      FluxInterpolatedBodyY[timeInterval][4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -0.5 * vzNodal[q];
      FluxInterpolatedBodyY[timeInterval][5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyY[timeInterval][6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -sxyNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyY[timeInterval][7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -syyNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyY[timeInterval][8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -syzNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyY[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      //--- z-dir
      FluxInterpolatedBodyZ[timeInterval][0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyZ[timeInterval][1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyZ[timeInterval][2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -vzNodal[q];
      FluxInterpolatedBodyZ[timeInterval][3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
      FluxInterpolatedBodyZ[timeInterval][4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -0.5 * vyNodal[q];
      FluxInterpolatedBodyZ[timeInterval][5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -0.5 * vxNodal[q];
      FluxInterpolatedBodyZ[timeInterval][6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -szxNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyZ[timeInterval][7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -syzNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyZ[timeInterval][8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
          -szzNodal[q] / data.material.local.rho;
      FluxInterpolatedBodyZ[timeInterval][9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0;
    }
  }
}

std::tuple<real, real, real, real, real, real, real, real, real, real, real, real, real>
    seissol::kernels::Time::calculateDamageAndBreakageStresses(real mu0,
                                                               real alpha,
                                                               real gammaR,
                                                               real xi0,
                                                               real xi,
                                                               real lambda0,
                                                               real EspI,
                                                               real EspII,
                                                               real stressXX,
                                                               real stressYY,
                                                               real stressZZ,
                                                               real stressXY,
                                                               real stressYZ,
                                                               real stressXZ,
                                                               real aB0,
                                                               real aB1,
                                                               real aB2,
                                                               real aB3) {
  real muEff = mu0 - alpha * gammaR * xi0 - 0.5 * alpha * gammaR * xi;
  real sxxS = lambda0 * EspI - alpha * gammaR * std::sqrt(EspII) + 2 * muEff * stressXX;
  real syyS = lambda0 * EspI - alpha * gammaR * std::sqrt(EspII) + 2 * muEff * stressYY;
  real szzS = lambda0 * EspI - alpha * gammaR * std::sqrt(EspII) + 2 * muEff * stressZZ;
  real sxyS = 2 * muEff * stressXY;
  real syzS = 2 * muEff * stressYZ;
  real szxS = 2 * muEff * stressXZ;
  real sxxB = (2.0 * aB2 + 3.0 * xi * aB3) * EspI + aB1 * std::sqrt(EspII) +
              (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressXX;
  real syyB = (2.0 * aB2 + 3.0 * xi * aB3) * EspI + aB1 * std::sqrt(EspII) +
              (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressYY;
  real szzB = (2.0 * aB2 + 3.0 * xi * aB3) * EspI + aB1 * std::sqrt(EspII) +
              (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressZZ;
  real sxyB = (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressXY;
  real syzB = (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressYZ;
  real szxB = (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressXZ;

  return std::make_tuple(
      muEff, sxxS, syyS, szzS, sxyS, syzS, szxS, sxxB, syyB, szzB, sxyB, syzB, szxB);
}

void seissol::kernels::Time::calculateStressesFromDamageAndBreakageStresses(real* qStressDofsXX,
                                                                            real* qStressDofsYY,
                                                                            real* qStressDofsZZ,
                                                                            real* qStressDofsXY,
                                                                            real* qStressDofsYZ,
                                                                            real* qStressDofsXZ,
                                                                            real* qStressDofsU,
                                                                            real* qStressDofsV,
                                                                            real* qStressDofsW,
                                                                            real* qStressDofsDAM,
                                                                            real* qStressDofsBRE,
                                                                            const real* qIU,
                                                                            const real* qIV,
                                                                            const real* qIW,
                                                                            const real* qIDAM,
                                                                            const real* qIBRE,
                                                                            real sxxS,
                                                                            real sxxB,
                                                                            real syyS,
                                                                            real syyB,
                                                                            real szzS,
                                                                            real szzB,
                                                                            real sxyS,
                                                                            real sxyB,
                                                                            real syzS,
                                                                            real syzB,
                                                                            real szxS,
                                                                            real szxB,
                                                                            unsigned int i) {
  qStressDofsXX[i] = (1 - qIBRE[i]) * sxxS + qIBRE[i] * sxxB;
  qStressDofsYY[i] = (1 - qIBRE[i]) * syyS + qIBRE[i] * syyB;
  qStressDofsZZ[i] = (1 - qIBRE[i]) * szzS + qIBRE[i] * szzB;
  qStressDofsXY[i] = (1 - qIBRE[i]) * sxyS + qIBRE[i] * sxyB;
  qStressDofsYZ[i] = (1 - qIBRE[i]) * syzS + qIBRE[i] * syzB;
  qStressDofsXZ[i] = (1 - qIBRE[i]) * szxS + qIBRE[i] * szxB;
  qStressDofsU[i] = qIU[i];
  qStressDofsV[i] = qIV[i];
  qStressDofsW[i] = qIW[i];
  qStressDofsDAM[i] = qIDAM[i];
  qStressDofsBRE[i] = qIBRE[i];
}

void seissol::kernels::Time::updateNonLinearMaterial(seissol::model::DamagedElasticMaterial& material){

}