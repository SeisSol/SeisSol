/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 * Dynamic Rupture kernel of SeisSol.
 **/

#include "DynamicRupture.h"

#ifndef NDEBUG
#pragma message "compiling dynamic rupture kernel with assertions"
#endif

#include <cassert>
#include <cstring>
#include <stdint.h>

#include <generated_code/kernel.h>
#include <Kernels/common.hpp>
#include <Numerical_aux/Quadrature.h>
#ifdef ACL_DEVICE
#include "device.h"
#endif
#include <yateto.h>

void seissol::kernels::DynamicRupture::checkGlobalData(GlobalData const* global, size_t alignment) {
#ifndef NDEBUG
  for (unsigned face = 0; face < 4; ++face) {
    for (unsigned h = 0; h < 4; ++h) {
      assert( ((uintptr_t const)global->faceToNodalMatrices(face, h)) % alignment == 0 );
    }
  }
#endif
}

void seissol::kernels::DynamicRupture::setHostGlobalData(GlobalData const* global) {
  checkGlobalData(global, ALIGNMENT);
  m_krnlPrototype.V3mTo2n = global->faceToNodalMatrices;
  m_timeKernel.setHostGlobalData(global);
}


void seissol::kernels::DynamicRupture::setGlobalData(const CompoundGlobalData& global) {
  this->setHostGlobalData(global.onHost);
#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  checkGlobalData(global.onDevice, deviceAlignment);
  m_gpuKrnlPrototype.V3mTo2n = global.onDevice->faceToNodalMatrices;
  m_timeKernel.setGlobalData(global);
#endif
}


void seissol::kernels::DynamicRupture::setTimeStepWidth(double timestep)
{
#ifdef USE_DR_CELLAVERAGE
  static_assert(false, "Cell average currently not supported");
  /*double subIntervalWidth = timestep / CONVERGENCE_ORDER;
  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    double t1 = timeInterval * subIntervalWidth;
    double t2 = t1 + subIntervalWidth;
    /// Compute time-integrated Taylor expansion (at t0=0) weights for interval [t1,t2].
    unsigned factorial = 1;
    for (unsigned derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
      m_timeFactors[timeInterval][derivative] = (t2-t1) / (factorial * subIntervalWidth);
      t1 *= t1;
      t2 *= t2;
      factorial *= (derivative+2);
    }
    /// We define the time "point" of the interval as the centre of the interval in order
    /// to be somewhat compatible to legacy code.
    timePoints[timeInterval] = timeInterval * subIntervalWidth + subIntervalWidth / 2.;
    timeWeights[timeInterval] = subIntervalWidth;
  }*/
#else
  seissol::quadrature::GaussLegendre(timePoints, timeWeights, CONVERGENCE_ORDER);
  for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
    timePoints[point] = 0.5 * (timestep * timePoints[point] + timestep);
    timeWeights[point] = 0.5 * timestep * timeWeights[point];
  }
#endif
}

void seissol::kernels::DynamicRupture::spaceTimeInterpolation(  DRFaceInformation const&    faceInfo,
                                                                GlobalData const*           global,
                                                                DRGodunovData const*        godunovData,
                                                                real const*                 timeDerivativePlus,
                                                                real const*                 timeDerivativeMinus,
                                                                real                        QInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()],
                                                                real                        QInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()],
                                                                real const*                 timeDerivativePlus_prefetch,
                                                                real const*                 timeDerivativeMinus_prefetch ) {
  // assert alignments
#ifndef NDEBUG
  assert( timeDerivativePlus != nullptr );
  assert( timeDerivativeMinus != nullptr );
  assert( ((uintptr_t)timeDerivativePlus) % ALIGNMENT == 0 );
  assert( ((uintptr_t)timeDerivativeMinus) % ALIGNMENT == 0 );
  assert( ((uintptr_t)&QInterpolatedPlus[0]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)&QInterpolatedMinus[0]) % ALIGNMENT == 0 );
  assert( tensor::Q::size() == tensor::I::size() );
#endif

  real degreesOfFreedomPlus[tensor::Q::size()] __attribute__((aligned(PAGESIZE_STACK)));
  real degreesOfFreedomMinus[tensor::Q::size()] __attribute__((aligned(PAGESIZE_STACK)));

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints krnl = m_krnlPrototype;

  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    m_timeKernel.computeTaylorExpansion(timePoints[timeInterval], 0.0, timeDerivativePlus, degreesOfFreedomPlus);
    m_timeKernel.computeTaylorExpansion(timePoints[timeInterval], 0.0, timeDerivativeMinus, degreesOfFreedomMinus);

    real const* plusPrefetch = (timeInterval < CONVERGENCE_ORDER-1) ? &QInterpolatedPlus[timeInterval+1][0] : timeDerivativePlus_prefetch;
    real const* minusPrefetch = (timeInterval < CONVERGENCE_ORDER-1) ? &QInterpolatedMinus[timeInterval+1][0] : timeDerivativeMinus_prefetch;
    
    krnl.QInterpolated = &QInterpolatedPlus[timeInterval][0];
    krnl.Q = degreesOfFreedomPlus;
    krnl.TinvT = godunovData->TinvT;
    krnl._prefetch.QInterpolated = plusPrefetch;
    krnl.execute(faceInfo.plusSide, 0);
    
    krnl.QInterpolated = &QInterpolatedMinus[timeInterval][0];
    krnl.Q = degreesOfFreedomMinus;
    krnl.TinvT = godunovData->TinvT;
    krnl._prefetch.QInterpolated = minusPrefetch;
    krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);
  }
}

void seissol::kernels::DynamicRupture::batchedSpaceTimeInterpolation(ConditionalBatchTableT& table) {
#ifdef ACL_DEVICE

  real** degreesOfFreedomPlus{nullptr};
  real** degreesOfFreedomMinus{nullptr};

  auto resetDeviceCurrentState = [this](size_t counter) {
    for (size_t i = 0; i < counter; ++i) {
      this->device.api->popStackMemory();
    }
    this->device.api->fastStreamsSync();
    this->device.api->resetCircularStreamCounter();
  };

  device.api->resetCircularStreamCounter();
  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    ConditionalKey timeIntegrationKey(*KernelNames::DrTime);
    if (table.find(timeIntegrationKey) != table.end()) {
      BatchTable &entry = table[timeIntegrationKey];

      unsigned maxNumElements = (entry.content[*EntityId::DrDerivativesPlus])->getSize();
      real** timeDerivativePlus = (entry.content[*EntityId::DrDerivativesPlus])->getPointers();
      degreesOfFreedomPlus = (entry.content[*EntityId::DrIdofsPlus])->getPointers();

      m_timeKernel.computeBatchedTaylorExpansion(timePoints[timeInterval],
                                                 0.0,
                                                 timeDerivativePlus,
                                                 degreesOfFreedomPlus,
                                                 maxNumElements);

      real** timeDerivativeMinus = (entry.content[*EntityId::DrDerivativesMinus])->getPointers();
      degreesOfFreedomMinus = (entry.content[*EntityId::DrIdofsMinus])->getPointers();
      m_timeKernel.computeBatchedTaylorExpansion(timePoints[timeInterval],
                                                 0.0,
                                                 timeDerivativeMinus,
                                                 degreesOfFreedomMinus,
                                                 maxNumElements);
    }

    size_t streamCounter{0};
    for (unsigned side = 0; side < 4; ++side) {
      ConditionalKey plusSideKey(*KernelNames::DrSpaceMap, side);
      if (table.find(plusSideKey) != table.end()) {
        BatchTable &entry = table[plusSideKey];
        const size_t numElements = (entry.content[*EntityId::DrIdofsPlus])->getSize();

        auto krnl = m_gpuKrnlPrototype;
        real *tmpMem = (real *) (device.api->getStackMemory(krnl.TmpMaxMemRequiredInBytes * numElements));
        ++streamCounter;
        krnl.linearAllocator.initialize(tmpMem);
        krnl.streamPtr = device.api->getNextCircularStream();
        krnl.numElements = numElements;

        krnl.QInterpolated = (entry.content[*EntityId::DrQInterpolatedPlus])->getPointers();
        krnl.extraOffset_QInterpolated = timeInterval * tensor::QInterpolated::size();
        krnl.Q = const_cast<real const **>((entry.content[*EntityId::DrIdofsPlus])->getPointers());
        krnl.TinvT = const_cast<real const **>((entry.content[*EntityId::DrTinvT])->getPointers());
        krnl.execute(side, 0);
      }

      for (unsigned faceRelation = 0; faceRelation < 4; ++faceRelation) {
        ConditionalKey minusSideKey(*KernelNames::DrSpaceMap, side, faceRelation);
        if (table.find(minusSideKey) != table.end()) {
          BatchTable &entry = table[minusSideKey];
          const size_t numElements = (entry.content[*EntityId::DrIdofsMinus])->getSize();

          auto krnl = m_gpuKrnlPrototype;
          real *tmpMem = (real *) (device.api->getStackMemory(krnl.TmpMaxMemRequiredInBytes * numElements));
          ++streamCounter;
          krnl.linearAllocator.initialize(tmpMem);
          krnl.streamPtr = device.api->getNextCircularStream();
          krnl.numElements = numElements;

          krnl.QInterpolated = (entry.content[*EntityId::DrQInterpolatedMinus])->getPointers();
          krnl.extraOffset_QInterpolated = timeInterval * tensor::QInterpolated::size();
          krnl.Q = const_cast<real const **>((entry.content[*EntityId::DrIdofsMinus])->getPointers());
          krnl.TinvT = const_cast<real const **>((entry.content[*EntityId::DrTinvT])->getPointers());
          krnl.execute(side, faceRelation);
        }
      }
    }
    resetDeviceCurrentState(streamCounter);
  }
  device.api->synchDevice();
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::DynamicRupture::flopsGodunovState( DRFaceInformation const&  faceInfo,
                                                          long long&                o_nonZeroFlops,
                                                          long long&                o_hardwareFlops )
{
  m_timeKernel.flopsTaylorExpansion(o_nonZeroFlops, o_hardwareFlops);
 
  // 2x evaluateTaylorExpansion
  o_nonZeroFlops *= 2;
  o_hardwareFlops *= 2;
  
  o_nonZeroFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::nonZeroFlops(faceInfo.plusSide, 0);
  o_hardwareFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::hardwareFlops(faceInfo.plusSide, 0);
  
  o_nonZeroFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::nonZeroFlops(faceInfo.minusSide, faceInfo.faceRelation);
  o_hardwareFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::hardwareFlops(faceInfo.minusSide, faceInfo.faceRelation);
  
  o_nonZeroFlops *= CONVERGENCE_ORDER;
  o_hardwareFlops *= CONVERGENCE_ORDER;
}
