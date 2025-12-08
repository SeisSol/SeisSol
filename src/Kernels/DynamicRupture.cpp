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
#include "Common/Marker.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Parallel/Runtime/Stream.h"

#include <Alignment.h>
#include <Initializer/Typedefs.h>
#include <Memory/GlobalData.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <cstring>
#include <stdint.h>
#include <utils/logger.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"

#include <Device/device.h>
#endif

#ifndef NDEBUG
#include <cstdint>
#endif

namespace seissol::kernels {

template <typename Cfg>
void DynamicRupture<Cfg>::setGlobalData(const GlobalData& global) {
  m_krnlPrototype.V3mTo2n = global.get<Cfg>().faceToNodalMatrices;
#ifdef ACL_DEVICE
  m_gpuKrnlPrototype.V3mTo2n = global.get<Cfg, Executor::Device>().faceToNodalMatrices;
#endif

  m_timeKernel.setGlobalData(global);
}

template <typename Cfg>
void DynamicRupture<Cfg>::spaceTimeInterpolation(
    const DRFaceInformation& faceInfo,
    const DRGodunovData<Cfg>* godunovData,
    const real* timeDerivativePlus,
    const real* timeDerivativeMinus,
    real qInterpolatedPlus[dr::misc::TimeSteps<Cfg>][seissol::tensor::QInterpolated<Cfg>::size()],
    real qInterpolatedMinus[dr::misc::TimeSteps<Cfg>][seissol::tensor::QInterpolated<Cfg>::size()],
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
  static_assert(tensor::Q<Cfg>::size() == tensor::I<Cfg>::size(),
                "The tensors Q and I need to match in size");
#endif

  alignas(PagesizeStack) real degreesOfFreedomPlus[tensor::Q<Cfg>::size()];
  alignas(PagesizeStack) real degreesOfFreedomMinus[tensor::Q<Cfg>::size()];

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg> krnl = m_krnlPrototype;
  for (std::size_t timeInterval = 0; timeInterval < dr::misc::TimeSteps<Cfg>; ++timeInterval) {
    m_timeKernel.evaluate(
        &coeffs[timeInterval * Cfg::ConvergenceOrder], timeDerivativePlus, degreesOfFreedomPlus);
    m_timeKernel.evaluate(
        &coeffs[timeInterval * Cfg::ConvergenceOrder], timeDerivativeMinus, degreesOfFreedomMinus);

    const real* plusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps<Cfg>)
                                   ? &qInterpolatedPlus[timeInterval + 1][0]
                                   : timeDerivativePlusPrefetch;
    const real* minusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps<Cfg>)
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

template <typename Cfg>
void DynamicRupture<Cfg>::batchedSpaceTimeInterpolation(
    SEISSOL_GPU_PARAM recording::DrConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  real** degreesOfFreedomPlus{nullptr};
  real** degreesOfFreedomMinus{nullptr};

  for (unsigned timeInterval = 0; timeInterval < dr::misc::TimeSteps<Cfg>; ++timeInterval) {
    ConditionalKey timeIntegrationKey(*KernelNames::DrTime);
    if (table.find(timeIntegrationKey) != table.end()) {
      auto& entry = table[timeIntegrationKey];

      unsigned maxNumElements = (entry.get(inner_keys::Dr::Id::DerivativesPlus))->getSize();
      real** timeDerivativePlus =
          (entry.get(inner_keys::Dr::Id::DerivativesPlus))->getDeviceDataPtrAs<real*>();
      degreesOfFreedomPlus =
          (entry.get(inner_keys::Dr::Id::IdofsPlus))->getDeviceDataPtrAs<real*>();

      m_timeKernel.evaluateBatched(&coeffs[timeInterval * Cfg::ConvergenceOrder],
                                   const_cast<const real**>(timeDerivativePlus),
                                   degreesOfFreedomPlus,
                                   maxNumElements,
                                   runtime);

      real** timeDerivativeMinus =
          (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtrAs<real*>();
      degreesOfFreedomMinus =
          (entry.get(inner_keys::Dr::Id::IdofsMinus))->getDeviceDataPtrAs<real*>();
      m_timeKernel.evaluateBatched(&coeffs[timeInterval * Cfg::ConvergenceOrder],
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
              (entry.get(inner_keys::Dr::Id::QInterpolatedPlus))->getDeviceDataPtrAs<real*>();
          krnl.extraOffset_QInterpolated = timeInterval * tensor::QInterpolated<Cfg>::size();
          krnl.Q = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::IdofsPlus))->getDeviceDataPtrAs<real*>());
          krnl.TinvT = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::TinvT))->getDeviceDataPtrAs<real*>());
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
              (entry.get(inner_keys::Dr::Id::QInterpolatedMinus))->getDeviceDataPtrAs<real*>();
          krnl.extraOffset_QInterpolated = timeInterval * tensor::QInterpolated<Cfg>::size();
          krnl.Q = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::IdofsMinus))->getDeviceDataPtrAs<real*>());
          krnl.TinvT = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::TinvT))->getDeviceDataPtrAs<real*>());
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

template <typename Cfg>
void DynamicRupture<Cfg>::flopsGodunovState(const DRFaceInformation& faceInfo,
                                            std::uint64_t& nonZeroFlops,
                                            std::uint64_t& hardwareFlops) {
  m_timeKernel.flopsEvaluate(nonZeroFlops, hardwareFlops);

  // 2x evaluateTaylorExpansion
  nonZeroFlops *= 2;
  hardwareFlops *= 2;

  nonZeroFlops +=
      dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg>::nonZeroFlops(
          faceInfo.plusSide, 0);
  hardwareFlops +=
      dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg>::hardwareFlops(
          faceInfo.plusSide, 0);

  nonZeroFlops +=
      dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg>::nonZeroFlops(
          faceInfo.minusSide, faceInfo.faceRelation);
  hardwareFlops +=
      dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg>::hardwareFlops(
          faceInfo.minusSide, faceInfo.faceRelation);

  nonZeroFlops *= dr::misc::TimeSteps<Cfg>;
  hardwareFlops *= dr::misc::TimeSteps<Cfg>;
}

#define SEISSOL_CONFIGITER(cfg) template class DynamicRupture<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels
