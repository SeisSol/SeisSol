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
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <cassert>
#include <cstring>
#include <iterator>
#include <stdint.h>
#include <yateto.h>
#include <yateto/InitTools.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"

#include <Device/device.h>
#endif

#ifndef ACL_DEVICE
#include <utils/logger.h>
#endif

#ifndef NDEBUG
#include <cstdint>
#endif

GENERATE_HAS_MEMBER(I)

namespace seissol::kernels {

#ifdef ACL_DEVICE
static_assert(*seissol::recording::DrFaceRelations::Count ==
              Cell::NumFaces * dr::misc::NumFaceRelations);
#endif

template <typename Cfg>
void DynamicRupture<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {
  krnlPrototype_.bindGlobals(*global.onHost);
#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  gpuKrnlPrototype_.bindGlobals(*global.onDevice);
  gpuCombinedKrnlPrototype_.bindGlobals(*global.onDevice);
#endif

  timeKernel_.setGlobalData(global);
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
  // The dynamic rupture families are indexed by the side and the face relation. Relation 0
  // addresses the plus side, relation 1 the minus side at a zero face orientation index, which the
  // canonical vertex numbering guarantees on every interior face.
  static_assert(std::size(dynamicRupture::kernel::nodalFlux<Cfg>::ExecutePtrs) ==
                Cell::NumFaces * dr::misc::NumFaceRelations);
  static_assert(
      std::size(
          dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg>::ExecutePtrs) ==
      Cell::NumFaces * dr::misc::NumFaceRelations);
  static_assert(std::size(tensor::V3mTo2n<Cfg>::Size) ==
                Cell::NumFaces * dr::misc::NumFaceRelations);
  static_assert(std::size(tensor::V3mTo2nTWDivM<Cfg>::Size) ==
                Cell::NumFaces * dr::misc::NumFaceRelations);

  // assert alignments
  assert(timeDerivativePlus != nullptr);
  assert(timeDerivativeMinus != nullptr);
  assert((reinterpret_cast<uintptr_t>(timeDerivativePlus)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivativeMinus)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedPlus[0])) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedMinus[0])) % Vectorsize == 0);
  static_assert(tensor::Q<Cfg>::size() == tensor::I<Cfg>::size(),
                "The tensors Q and I need to match in size");

  alignas(PagesizeStack) real degreesOfFreedomPlus[tensor::Q<Cfg>::size()];
  alignas(PagesizeStack) real degreesOfFreedomMinus[tensor::Q<Cfg>::size()];

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg> krnl = krnlPrototype_;
  for (std::size_t timeInterval = 0; timeInterval < dr::misc::TimeSteps<Cfg>; ++timeInterval) {
    timeKernel_.evaluate(
        &coeffs[timeInterval * Cfg::ConvergenceOrder], timeDerivativePlus, degreesOfFreedomPlus);
    timeKernel_.evaluate(
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

  // interpolate all timesteps in a single kernel

  runtime.envMany(Cell::NumFaces * dr::misc::NumFaceRelations, [&](void* stream, size_t i) {
    const auto side = i / dr::misc::NumFaceRelations;
    const auto faceRelation = i % dr::misc::NumFaceRelations;

    const ConditionalKey minusSideKey(*KernelNames::DrSpaceMap, side, faceRelation);
    if (table.find(minusSideKey) != table.end()) {
      auto& entry = table[minusSideKey];
      const size_t numElements = (entry.get<real*>(inner_keys::Dr::Id::IdofsMinus))->getSize();

      auto krnl = gpuCombinedKrnlPrototype_;
      real* tmpMem = reinterpret_cast<real*>(
          device_.api().allocMemAsync(krnl.TmpMaxMemRequiredInBytes * numElements, stream));
      krnl.linearAllocator.initialize(tmpMem);
      krnl.streamPtr = stream;
      krnl.numElements = numElements;

      std::size_t offsetQDR = 0;
      for (std::size_t s = 0; s < dr::misc::TimeSteps<Cfg>; ++s) {
        krnl.QDR(s) =
            (entry.get<real*>(inner_keys::Dr::Id::QInterpolatedMinus))->getDeviceDataPtr();
        krnl.extraOffset_QDR(s) = offsetQDR;
        offsetQDR += tensor::QDR<Cfg>::size(s);
      }

      std::size_t offsetDQ = 0;
      for (std::size_t p = 0; p < Cfg::ConvergenceOrder; ++p) {
        krnl.dQ(p) = const_cast<const real**>(
            (entry.get<real*>(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr());
        krnl.extraOffset_dQ(p) = offsetDQ;
        offsetDQ += tensor::dQ<Cfg>::size(p);
      }

      for (std::size_t s = 0; s < dr::misc::TimeSteps<Cfg>; ++s) {
        for (std::size_t p = 0; p < Cfg::ConvergenceOrder; ++p) {
          krnl.coeffDR(s * Cfg::ConvergenceOrder + p) = coeffs[s * Cfg::ConvergenceOrder + p];
        }
      }

      set_I(krnl, entry.get<real*>(inner_keys::Dr::Id::IdofsMinus)->getDeviceDataPtr());

      krnl.TinvT = const_cast<const real**>(
          (entry.get<real*>(inner_keys::Dr::Id::TinvT))->getDeviceDataPtr());
      krnl.execute(side, faceRelation);

      device_.api().freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
    }
  });
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
PerformanceEstimate DynamicRupture<Cfg>::metrics(const DRFaceInformation& faceInfo) const {
  if (isDeviceOn()) {
    return PerformanceEstimate::fromKernel<dynamicRupture::kernel::projectToDR<Cfg>>(
               faceInfo.plusSide, 0) +
           PerformanceEstimate::fromKernel<dynamicRupture::kernel::projectToDR<Cfg>>(
               faceInfo.minusSide, faceInfo.faceRelation);
  } else {
    auto estimate = timeKernel_.metrics();

    // 2x evaluateTaylorExpansion
    estimate *= 2;

    estimate += PerformanceEstimate::fromKernel<
        dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg>>(faceInfo.plusSide, 0);

    estimate += PerformanceEstimate::fromKernel<
        dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg>>(
        faceInfo.minusSide, faceInfo.faceRelation);

    estimate *= dr::misc::TimeSteps<Cfg>;

    // legacy CPU memory estimate
    estimate.bytes = (tensor::TinvT<Cfg>::size() +
                      tensor::QInterpolated<Cfg>::size() * 2 * dr::misc::TimeSteps<Cfg> +
                      yateto::computeFamilySize<tensor::dQ<Cfg>>() * 2) *
                     sizeof(real);

    return estimate;
  }
}

#define SEISSOL_INSTANTIATE(Cfg) template class DynamicRupture<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
