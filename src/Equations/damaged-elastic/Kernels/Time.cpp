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
    std::tie(EspI, EspII, xi) = calculateEsp(
        exxNodal, eyyNodal, ezzNodal, exyNodal, eyzNodal, ezxNodal, q, m_damagedElasticParameters);

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
  real xi = computexi(EspI, EspII);
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

    real mu_eff;
    Stresses sSp, sBp;

    std::tie(mu_eff, sSp, sBp) = calculateDamageAndBreakageStresses(
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
    seissol::kernels::Time::StressDofs sDofsPlus;
    sDofsPlus.qStressDofsXX = &dofsStressNPlus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsYY = &dofsStressNPlus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsZZ = &dofsStressNPlus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsXY = &dofsStressNPlus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsYZ = &dofsStressNPlus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsXZ = &dofsStressNPlus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsU = &dofsStressNPlus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsV = &dofsStressNPlus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsW = &dofsStressNPlus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsDAM = &dofsStressNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsPlus.qStressDofsBRE = &dofsStressNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];

    calculateStressesFromDamageAndBreakageStresses(
        sDofsPlus,
        &dofsNPlus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNPlus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNPlus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNPlus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNPlus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        sSp,
        sBp,
        q);

    std::tie(EspIm, EspIIm, xim) = calculateEsp(&dofsNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                                                &dofsNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                                                &dofsNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                                                &dofsNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                                                &dofsNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                                                &dofsNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                                                q,
                                                &damagedElasticParameters);
    real alpham = dofsNMinus[9];

    Stresses sSm, sBm;
    std::tie(mu_eff, sSm, sBm) = calculateDamageAndBreakageStresses(
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

    seissol::kernels::Time::StressDofs sDofsMinus;
    sDofsMinus.qStressDofsXX = &dofsStressNMinus[0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsYY = &dofsStressNMinus[1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsZZ = &dofsStressNMinus[2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsXY = &dofsStressNMinus[3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsYZ = &dofsStressNMinus[4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsXZ = &dofsStressNMinus[5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsU = &dofsStressNMinus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsV = &dofsStressNMinus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsW = &dofsStressNMinus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsDAM = &dofsStressNMinus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    sDofsMinus.qStressDofsBRE = &dofsStressNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];

    calculateStressesFromDamageAndBreakageStresses(
        sDofsMinus,
        &dofsNMinus[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNMinus[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNMinus[8 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNMinus[9 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        &dofsNMinus[10 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
        sSm,
        sBm,
        q);
  }
}

void seissol::kernels::Time::stressToDofsDynamicRupture(real* dofsStress, const real* dofs) {
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

  const real aB0 = m_damagedElasticParameters->aB0;
  const real aB1 = m_damagedElasticParameters->aB1;
  const real aB2 = m_damagedElasticParameters->aB2;
  const real aB3 = m_damagedElasticParameters->aB3;
  auto getQ = [](const real* qI, unsigned o, unsigned q) {
    constexpr size_t offset1 =
        seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints;
    constexpr size_t offset2 = seissol::dr::misc::numPaddedPoints;
    return &qI[o * offset1 + q * offset2];
  };

  auto getQStress = [](real* qStressI, unsigned o, unsigned q) {
    constexpr size_t offset1 =
        seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints;
    constexpr size_t offset2 = seissol::dr::misc::numPaddedPoints;
    return &qStressI[o * offset1 + q * offset2];
  };

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints; i++) {

      const real EspIp =
          (getQ(qIPlus, o, XX)[i]) + (getQ(qIPlus, o, YY)[i]) + (getQ(qIPlus, o, ZZ)[i]);
      const real EspIIp = (getQ(qIPlus, o, XX)[i]) * (getQ(qIPlus, o, XX)[i]) +
                          (getQ(qIPlus, o, YY)[i]) * (getQ(qIPlus, o, YY)[i]) +
                          (getQ(qIPlus, o, ZZ)[i]) * (getQ(qIPlus, o, ZZ)[i]) +
                          2 * (getQ(qIPlus, o, XY)[i]) * (getQ(qIPlus, o, XY)[i]) +
                          2 * (getQ(qIPlus, o, YZ)[i]) * (getQ(qIPlus, o, YZ)[i]) +
                          2 * (getQ(qIPlus, o, XZ)[i]) * (getQ(qIPlus, o, XZ)[i]);
      real alphap = getQ(qIPlus, o, DAM)[i];
      real xip = computexi(EspIp, EspIIp);

      real mu_eff;
      Stresses sSp, sBp;
      std::tie(mu_eff, sSp, sBp) = calculateDamageAndBreakageStresses(mu0P,
                                                                      alphap,
                                                                      impAndEta[ltsFace].gammaRP,
                                                                      impAndEta[ltsFace].xi0P,
                                                                      xip,
                                                                      lambda0P,
                                                                      EspIp,
                                                                      EspIIp,
                                                                      getQ(qIPlus, o, XX)[i],
                                                                      getQ(qIPlus, o, YY)[i],
                                                                      getQ(qIPlus, o, ZZ)[i],
                                                                      getQ(qIPlus, o, XY)[i],
                                                                      getQ(qIPlus, o, YZ)[i],
                                                                      getQ(qIPlus, o, XZ)[i],
                                                                      aB0,
                                                                      aB1,
                                                                      aB2,
                                                                      aB3);
      seissol::kernels::Time::StressDofs sDofsPlus;
      sDofsPlus.qStressDofsXX = getQStress(qStressIPlus, o, XX);
      sDofsPlus.qStressDofsYY = getQStress(qStressIPlus, o, YY);
      sDofsPlus.qStressDofsZZ = getQStress(qStressIPlus, o, ZZ);
      sDofsPlus.qStressDofsXY = getQStress(qStressIPlus, o, XY);
      sDofsPlus.qStressDofsYZ = getQStress(qStressIPlus, o, YZ);
      sDofsPlus.qStressDofsXZ = getQStress(qStressIPlus, o, XZ);
      sDofsPlus.qStressDofsU = getQStress(qStressIPlus, o, U);
      sDofsPlus.qStressDofsV = getQStress(qStressIPlus, o, V);
      sDofsPlus.qStressDofsW = getQStress(qStressIPlus, o, W);
      sDofsPlus.qStressDofsDAM = getQStress(qStressIPlus, o, DAM);
      sDofsPlus.qStressDofsBRE = getQStress(qStressIPlus, o, BRE);

      calculateStressesFromDamageAndBreakageStresses(sDofsPlus,
                                                     getQ(qIPlus, o, U),
                                                     getQ(qIPlus, o, V),
                                                     getQ(qIPlus, o, W),
                                                     getQ(qIPlus, o, DAM),
                                                     getQ(qIPlus, o, BRE),
                                                     sSp,
                                                     sBp,
                                                     i);

      const real EspIm =
          (getQ(qIMinus, o, XX)[i]) + (getQ(qIMinus, o, YY)[i]) + (getQ(qIMinus, o, ZZ)[i]);
      const real EspIIm = (getQ(qIMinus, o, XX)[i]) * (getQ(qIMinus, o, XX)[i]) +
                          (getQ(qIMinus, o, YY)[i]) * (getQ(qIMinus, o, YY)[i]) +
                          (getQ(qIMinus, o, ZZ)[i]) * (getQ(qIMinus, o, ZZ)[i]) +
                          2 * (getQ(qIMinus, o, XY)[i]) * (getQ(qIMinus, o, XY)[i]) +
                          2 * (getQ(qIMinus, o, YZ)[i]) * (getQ(qIMinus, o, YZ)[i]) +
                          2 * (getQ(qIMinus, o, XZ)[i]) * (getQ(qIMinus, o, XZ)[i]);

      real alpham = getQ(qIMinus, o, DAM)[i];
      real xim = computexi(EspIm, EspIIm);

      Stresses sSm, sBm;

      std::tie(mu_eff, sSm, sBm) = calculateDamageAndBreakageStresses(mu0M,
                                                                      alpham,
                                                                      impAndEta[ltsFace].gammaRM,
                                                                      impAndEta[ltsFace].xi0M,
                                                                      xim,
                                                                      lambda0M,
                                                                      EspIm,
                                                                      EspIIm,
                                                                      getQ(qIMinus, o, XX)[i],
                                                                      getQ(qIMinus, o, YY)[i],
                                                                      getQ(qIMinus, o, ZZ)[i],
                                                                      getQ(qIMinus, o, XY)[i],
                                                                      getQ(qIMinus, o, YZ)[i],
                                                                      getQ(qIMinus, o, XZ)[i],
                                                                      aB0,
                                                                      aB1,
                                                                      aB2,
                                                                      aB3);

      seissol::kernels::Time::StressDofs sDofsMinus;
      sDofsMinus.qStressDofsXX = getQStress(qStressIMinus, o, XX);
      sDofsMinus.qStressDofsYY = getQStress(qStressIMinus, o, YY);
      sDofsMinus.qStressDofsZZ = getQStress(qStressIMinus, o, ZZ);
      sDofsMinus.qStressDofsXY = getQStress(qStressIMinus, o, XY);
      sDofsMinus.qStressDofsYZ = getQStress(qStressIMinus, o, YZ);
      sDofsMinus.qStressDofsXZ = getQStress(qStressIMinus, o, XZ);
      sDofsMinus.qStressDofsU = getQStress(qStressIMinus, o, U);
      sDofsMinus.qStressDofsV = getQStress(qStressIMinus, o, V);
      sDofsMinus.qStressDofsW = getQStress(qStressIMinus, o, W);
      sDofsMinus.qStressDofsDAM = getQStress(qStressIMinus, o, DAM);
      sDofsMinus.qStressDofsBRE = getQStress(qStressIMinus, o, BRE);

      calculateStressesFromDamageAndBreakageStresses(sDofsMinus,
                                                     getQ(qIMinus, o, U),
                                                     getQ(qIMinus, o, V),
                                                     getQ(qIMinus, o, W),
                                                     getQ(qIMinus, o, DAM),
                                                     getQ(qIMinus, o, BRE),
                                                     sSm,
                                                     sBm,
                                                     i);
    }
  } // time integration loop
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
      real xi = computexi(EspI, EspII);

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

      real mu_eff;
      Stresses sS, sB;
      std::tie(mu_eff, sS, sB) = calculateDamageAndBreakageStresses(data.material.local.mu0,
                                                                    alphaNodal[q],
                                                                    data.material.local.gammaR,
                                                                    data.material.local.xi0,
                                                                    xi,
                                                                    data.material.local.lambda0,
                                                                    EspI,
                                                                    EspII,
                                                                    exxNodal[q] + epsInitxx,
                                                                    eyyNodal[q] + epsInityy,
                                                                    ezzNodal[q] + epsInitzz,
                                                                    exyNodal[q] + epsInitxy,
                                                                    eyzNodal[q] + epsInityz,
                                                                    ezxNodal[q] + epsInitzx,
                                                                    aB0,
                                                                    aB1,
                                                                    aB2,
                                                                    aB3);

      sxxNodal[q] = (1 - breakNodal[q]) * sS.sxx + breakNodal[q] * sB.sxx;
      syyNodal[q] = (1 - breakNodal[q]) * sS.syy + breakNodal[q] * sS.syy;
      szzNodal[q] = (1 - breakNodal[q]) * sS.szz + breakNodal[q] * sS.szz;
      sxyNodal[q] = (1 - breakNodal[q]) * sS.sxy + breakNodal[q] * sS.sxy;
      syzNodal[q] = (1 - breakNodal[q]) * sS.syz + breakNodal[q] * sS.syz;
      szxNodal[q] = (1 - breakNodal[q]) * sS.sxz + breakNodal[q] * sS.sxz;

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

std::tuple<real, seissol::kernels::Time::Stresses, seissol::kernels::Time::Stresses>
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
  seissol::kernels::Time::Stresses sS, sB;
  real muEff = mu0 - alpha * gammaR * xi0 - 0.5 * alpha * gammaR * xi;
  sS.sxx = lambda0 * EspI - alpha * gammaR * std::sqrt(EspII) + 2 * muEff * stressXX;
  sS.syy = lambda0 * EspI - alpha * gammaR * std::sqrt(EspII) + 2 * muEff * stressYY;
  sS.szz = lambda0 * EspI - alpha * gammaR * std::sqrt(EspII) + 2 * muEff * stressZZ;
  sS.sxy = 2 * muEff * stressXY;
  sS.syz = 2 * muEff * stressYZ;
  sS.sxz = 2 * muEff * stressXZ;
  sB.sxx = (2.0 * aB2 + 3.0 * xi * aB3) * EspI + aB1 * std::sqrt(EspII) +
           (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressXX;
  sB.syy = (2.0 * aB2 + 3.0 * xi * aB3) * EspI + aB1 * std::sqrt(EspII) +
           (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressYY;
  sB.szz = (2.0 * aB2 + 3.0 * xi * aB3) * EspI + aB1 * std::sqrt(EspII) +
           (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressZZ;
  sB.sxy = (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressXY;
  sB.syz = (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressYZ;
  sB.sxz = (2.0 * aB0 + aB1 * xi - aB3 * xi * xi * xi) * stressXZ;

  return std::make_tuple(muEff, sS, sB);
}

void seissol::kernels::Time::calculateStressesFromDamageAndBreakageStresses(
    seissol::kernels::Time::StressDofs& stressDofs,
    const real* qIU,
    const real* qIV,
    const real* qIW,
    const real* qIDAM,
    const real* qIBRE,
    const Stresses& sS,
    const Stresses& sB,
    unsigned int i) {
  stressDofs.qStressDofsXX[i] = (1 - qIBRE[i]) * sS.sxx + qIBRE[i] * sB.sxx;
  stressDofs.qStressDofsYY[i] = (1 - qIBRE[i]) * sS.syy + qIBRE[i] * sB.syy;
  stressDofs.qStressDofsZZ[i] = (1 - qIBRE[i]) * sS.szz + qIBRE[i] * sB.szz;
  stressDofs.qStressDofsXY[i] = (1 - qIBRE[i]) * sS.sxy + qIBRE[i] * sB.sxy;
  stressDofs.qStressDofsYZ[i] = (1 - qIBRE[i]) * sS.syz + qIBRE[i] * sB.syz;
  stressDofs.qStressDofsXZ[i] = (1 - qIBRE[i]) * sS.sxz + qIBRE[i] * sS.sxz;
  stressDofs.qStressDofsU[i] = qIU[i];
  stressDofs.qStressDofsV[i] = qIV[i];
  stressDofs.qStressDofsW[i] = qIW[i];
  stressDofs.qStressDofsDAM[i] = qIDAM[i];
  stressDofs.qStressDofsBRE[i] = qIBRE[i];
}

void seissol::kernels::Time::updateNonLinearMaterial(
    seissol::model::DamagedElasticMaterial& material, const real* Q_aveData) {
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

  real xi = computexi(EspI, EspII);

  const real alphaAve = Q_aveData[9];
  const real breakAve = Q_aveData[10];

  real lambda0 = material.lambda0;
  real mu0 = material.mu0;
  const real aB0 = m_damagedElasticParameters->aB0;
  const real aB1 = m_damagedElasticParameters->aB1;
  const real aB2 = m_damagedElasticParameters->aB2;
  const real aB3 = m_damagedElasticParameters->aB3;
  material.mu = (1 - breakAve) * (mu0 - alphaAve * material.xi0 * material.gammaR -
                                  0.5 * alphaAve * material.gammaR * xi) +
                breakAve * ((aB0 + 0.5 * aB1 * xi - 0.5 * aB3 * xi * xi * xi));
  material.lambda =
      (1 - breakAve) *
          (lambda0 - alphaAve * material.gammaR * (Q_aveData[0] + epsInitxx) / std::sqrt(EspII)) +
      breakAve *
          ((2.0 * aB2 + 3.0 * aB3 * xi) + aB1 * (Q_aveData[0] + epsInitxx) / std::sqrt(EspII));
  material.gamma = alphaAve * material.gammaR;
  material.epsxx_alpha = (Q_aveData[0] + epsInitxx);
  material.epsyy_alpha = (Q_aveData[1] + epsInityy);
  material.epszz_alpha = (Q_aveData[2] + epsInitzz);
  material.epsxy_alpha = (Q_aveData[3] + epsInitxy);
  material.epsyz_alpha = (Q_aveData[4] + epsInityz);
  material.epszx_alpha = (Q_aveData[5] + epsInitzx);
}

real seissol::kernels::Time::computexi(real EspI, real EspII) {
  return (EspII > 1e-30) ? EspI / std::sqrt(EspII) : 0.0;
}

std::tuple<real, real, real> seissol::kernels::Time::computealphalambdamu(const real* q,
                                                                          unsigned int o,
                                                                          unsigned int i,
                                                                          real lambda0,
                                                                          real mu0,
                                                                          real gammaR,
                                                                          real epsInitxx,
                                                                          real EspII,
                                                                          real aB0,
                                                                          real aB1,
                                                                          real aB2,
                                                                          real aB3,
                                                                          real xi,
                                                                          real xi0) {
  using namespace seissol::dr::misc::quantity_indices;

  auto getQ = [](const real* qI, unsigned o, unsigned q) {
    constexpr size_t offset1 =
        seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints;
    constexpr size_t offset2 = seissol::dr::misc::numPaddedPoints;
    return &qI[o * offset1 + q * offset2];
  };
  real alpha = getQ(q, o, DAM)[i];
  real lambda =
      (1 - getQ(q, o, BRE)[i]) *
          (lambda0 - alpha * gammaR * (getQ(q, o, XX)[i] + epsInitxx) / std::sqrt(EspII)) +
      getQ(q, o, BRE)[i] *
          (2.0 * aB2 + 3.0 * xi * aB3 + aB1 * (getQ(q, o, XX)[i] + epsInitxx) / std::sqrt(EspII));
  real mu = (1 - getQ(q, o, BRE)[i]) * (mu0 - alpha * xi0 * gammaR - 0.5 * alpha * gammaR * xi) +
            getQ(q, o, BRE)[i] * (aB0 + 0.5 * xi * aB1 - 0.5 * xi * xi * xi * aB3);

  return std::make_tuple(alpha, lambda, mu);
}
