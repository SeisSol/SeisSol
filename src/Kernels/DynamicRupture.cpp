// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "DynamicRupture.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "DataTypes/ConditionalTable.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Parallel/Runtime/Stream.h"
#include "utils/logger.h"

#include <cassert>
#include <cstring>
#include <stdint.h>

#ifdef ACL_DEVICE
#include "DataTypes/ConditionalKey.h"
#include "DataTypes/EncodedConstants.h"
#include "device.h"
#endif
#include <yateto.h>

#ifndef NDEBUG
#include <cstdint>
#endif

namespace seissol::kernels {

void DynamicRupture::setGlobalData(const CompoundGlobalData& global) {
  m_krnlPrototype.V3mTo2n = global.onHost->faceToNodalMatrices;
#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  m_gpuKrnlPrototype.V3mTo2n = global.onDevice->faceToNodalMatrices;
#endif

  m_timeKernel.setGlobalData(global);
}

void DynamicRupture::spaceTimeInterpolation(
    const DRFaceInformation& faceInfo,
    const GlobalData* global,
    const DRGodunovData* godunovData,
    DREnergyOutput* drEnergyOutput,
    const real* timeDerivativePlus,
    const real* timeDerivativeMinus,
    real qInterpolatedPlus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
    real qInterpolatedMinus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
    const real* timeDerivativePlusPrefetch,
    const real* timeDerivativeMinusPrefetch,
    const real* coeffs) {
  // assert alignments
#ifndef NDEBUG
  assert(timeDerivativePlus != nullptr);
  assert(timeDerivativeMinus != nullptr);
  assert((reinterpret_cast<uintptr_t>(timeDerivativePlus)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivativeMinus)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedPlus[0])) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedMinus[0])) % Alignment == 0);
  static_assert(tensor::Q::size() == tensor::I::size(),
                "The tensors Q and I need to match in size");
#endif

  alignas(PagesizeStack) real degreesOfFreedomPlus[tensor::Q::size()];
  alignas(PagesizeStack) real degreesOfFreedomMinus[tensor::Q::size()];

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints krnl = m_krnlPrototype;
  for (std::size_t timeInterval = 0; timeInterval < dr::misc::TimeSteps; ++timeInterval) {
    m_timeKernel.evaluate(
        &coeffs[timeInterval * ConvergenceOrder], timeDerivativePlus, degreesOfFreedomPlus);
    m_timeKernel.evaluate(
        &coeffs[timeInterval * ConvergenceOrder], timeDerivativeMinus, degreesOfFreedomMinus);

    const real* plusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps)
                                   ? &qInterpolatedPlus[timeInterval + 1][0]
                                   : timeDerivativePlusPrefetch;
    const real* minusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps)
                                    ? &qInterpolatedMinus[timeInterval + 1][0]
                                    : timeDerivativeMinusPrefetch;

    krnl.QInterpolated = &qInterpolatedPlus[timeInterval][0];
    krnl.Q = degreesOfFreedomPlus;
    krnl.TinvT = godunovData->dataTinvT;
    krnl._prefetch.QInterpolated = plusPrefetch;
    krnl.execute(faceInfo.plusSide, 0);

    krnl.QInterpolated = &qInterpolatedMinus[timeInterval][0];
    krnl.Q = degreesOfFreedomMinus;
    krnl.TinvT = godunovData->dataTinvT;
    krnl._prefetch.QInterpolated = minusPrefetch;
    krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);
  }
}

void DynamicRupture::batchedSpaceTimeInterpolation(
    DrConditionalPointersToRealsTable& table,
    const real* coeffs,
    seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  real** degreesOfFreedomPlus{nullptr};
  real** degreesOfFreedomMinus{nullptr};

  for (unsigned timeInterval = 0; timeInterval < dr::misc::TimeSteps; ++timeInterval) {
    ConditionalKey timeIntegrationKey(*KernelNames::DrTime);
    if (table.find(timeIntegrationKey) != table.end()) {
      auto& entry = table[timeIntegrationKey];

      unsigned maxNumElements = (entry.get(inner_keys::Dr::Id::DerivativesPlus))->getSize();
      real** timeDerivativePlus =
          (entry.get(inner_keys::Dr::Id::DerivativesPlus))->getDeviceDataPtr();
      degreesOfFreedomPlus = (entry.get(inner_keys::Dr::Id::IdofsPlus))->getDeviceDataPtr();

      m_timeKernel.evaluateBatched(&coeffs[timeInterval * ConvergenceOrder],
                                   const_cast<const real**>(timeDerivativePlus),
                                   degreesOfFreedomPlus,
                                   maxNumElements,
                                   runtime);

      real** timeDerivativeMinus =
          (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr();
      degreesOfFreedomMinus = (entry.get(inner_keys::Dr::Id::IdofsMinus))->getDeviceDataPtr();
      m_timeKernel.evaluateBatched(&coeffs[timeInterval * ConvergenceOrder],
                                   const_cast<const real**>(timeDerivativeMinus),
                                   degreesOfFreedomMinus,
                                   maxNumElements,
                                   runtime);
    }

    // finish all previous work in the default stream
    size_t streamCounter{0};
    runtime.envMany(20, [&](void* stream, size_t i) {
      unsigned side = i / 5;
      unsigned faceRelation = i % 5;
      if (faceRelation == 4) {
        ConditionalKey plusSideKey(*KernelNames::DrSpaceMap, side);
        if (table.find(plusSideKey) != table.end()) {
          auto& entry = table[plusSideKey];
          const size_t numElements = (entry.get(inner_keys::Dr::Id::IdofsPlus))->getSize();

          auto krnl = m_gpuKrnlPrototype;
          real* tmpMem = reinterpret_cast<real*>(
              device.api->allocMemAsync(krnl.TmpMaxMemRequiredInBytes * numElements, stream));
          ++streamCounter;
          krnl.linearAllocator.initialize(tmpMem);
          krnl.streamPtr = stream;
          krnl.numElements = numElements;

          krnl.QInterpolated =
              (entry.get(inner_keys::Dr::Id::QInterpolatedPlus))->getDeviceDataPtr();
          krnl.extraOffset_QInterpolated = timeInterval * tensor::QInterpolated::size();
          krnl.Q = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::IdofsPlus))->getDeviceDataPtr());
          krnl.TinvT =
              const_cast<const real**>((entry.get(inner_keys::Dr::Id::TinvT))->getDeviceDataPtr());
          krnl.execute(side, 0);

          device.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
        }
      } else {
        ConditionalKey minusSideKey(*KernelNames::DrSpaceMap, side, faceRelation);
        if (table.find(minusSideKey) != table.end()) {
          auto& entry = table[minusSideKey];
          const size_t numElements = (entry.get(inner_keys::Dr::Id::IdofsMinus))->getSize();

          auto krnl = m_gpuKrnlPrototype;
          real* tmpMem = reinterpret_cast<real*>(
              device.api->allocMemAsync(krnl.TmpMaxMemRequiredInBytes * numElements, stream));
          ++streamCounter;
          krnl.linearAllocator.initialize(tmpMem);
          krnl.streamPtr = stream;
          krnl.numElements = numElements;

          krnl.QInterpolated =
              (entry.get(inner_keys::Dr::Id::QInterpolatedMinus))->getDeviceDataPtr();
          krnl.extraOffset_QInterpolated = timeInterval * tensor::QInterpolated::size();
          krnl.Q = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::IdofsMinus))->getDeviceDataPtr());
          krnl.TinvT =
              const_cast<const real**>((entry.get(inner_keys::Dr::Id::TinvT))->getDeviceDataPtr());
          krnl.execute(side, faceRelation);

          device.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
        }
      }
    });
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void DynamicRupture::flopsGodunovState(const DRFaceInformation& faceInfo,
                                       std::uint64_t& nonZeroFlops,
                                       std::uint64_t& hardwareFlops) {
  m_timeKernel.flopsEvaluate(nonZeroFlops, hardwareFlops);

  // 2x evaluateTaylorExpansion
  nonZeroFlops *= 2;
  hardwareFlops *= 2;

  nonZeroFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::nonZeroFlops(
      faceInfo.plusSide, 0);
  hardwareFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::hardwareFlops(
      faceInfo.plusSide, 0);

  nonZeroFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::nonZeroFlops(
      faceInfo.minusSide, faceInfo.faceRelation);
  hardwareFlops += dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints::hardwareFlops(
      faceInfo.minusSide, faceInfo.faceRelation);

  nonZeroFlops *= dr::misc::TimeSteps;
  hardwareFlops *= dr::misc::TimeSteps;
}

} // namespace seissol::kernels
