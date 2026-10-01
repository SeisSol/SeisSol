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
#include "Kernels/Precision.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <cassert>
#include <cstring>
#include <iterator>
#include <stdint.h>
#include <utils/logger.h>
#include <yateto.h>
#include <yateto/InitTools.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"

#include <Device/device.h>
#endif

#ifndef NDEBUG
#include <cstdint>
#endif

GENERATE_HAS_MEMBER(I)

namespace seissol::kernels {

// The dynamic rupture families are indexed by the side and the face relation. Relation 0
// addresses the plus side, relation 1 the minus side at a zero face orientation index, which the
// canonical vertex numbering guarantees on every interior face.
static_assert(std::size(dynamicRupture::kernel::nodalFlux<Config>::ExecutePtrs) ==
              Cell::NumFaces * dr::misc::NumFaceRelations);
static_assert(
    std::size(
        dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Config>::ExecutePtrs) ==
    Cell::NumFaces * dr::misc::NumFaceRelations);
static_assert(std::size(tensor::V3mTo2n<Config>::Size) ==
              Cell::NumFaces * dr::misc::NumFaceRelations);
static_assert(std::size(tensor::V3mTo2nTWDivM<Config>::Size) ==
              Cell::NumFaces * dr::misc::NumFaceRelations);

#ifdef ACL_DEVICE
static_assert(*seissol::recording::DrFaceRelations::Count ==
              Cell::NumFaces * dr::misc::NumFaceRelations);
#endif

void DynamicRupture::setGlobalData(const CompoundGlobalData& global) {
  krnlPrototype_.bindGlobals(*global.onHost);
#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  gpuKrnlPrototype_.bindGlobals(*global.onDevice);
  gpuCombinedKrnlPrototype_.bindGlobals(*global.onDevice);
#endif

  timeKernel_.setGlobalData(global);
}

void DynamicRupture::spaceTimeInterpolation(
    const DRFaceInformation& faceInfo,
    const DRGodunovData* godunovData,
    const real* timeDerivativePlus,
    const real* timeDerivativeMinus,
    real qInterpolatedPlus[dr::misc::TimeSteps<Config>]
                          [seissol::tensor::QInterpolated<Config>::size()],
    real qInterpolatedMinus[dr::misc::TimeSteps<Config>]
                           [seissol::tensor::QInterpolated<Config>::size()],
    const real* timeDerivativePlusPrefetch,
    const real* timeDerivativeMinusPrefetch,
    const real* coeffs) {

  // assert alignments
  assert(timeDerivativePlus != nullptr);
  assert(timeDerivativeMinus != nullptr);
  assert((reinterpret_cast<uintptr_t>(timeDerivativePlus)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivativeMinus)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedPlus[0])) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedMinus[0])) % Vectorsize == 0);
  static_assert(tensor::Q<Config>::size() == tensor::I<Config>::size(),
                "The tensors Q and I need to match in size");

  alignas(PagesizeStack) real degreesOfFreedomPlus[tensor::Q<Config>::size()];
  alignas(PagesizeStack) real degreesOfFreedomMinus[tensor::Q<Config>::size()];

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Config> krnl = krnlPrototype_;
  for (std::size_t timeInterval = 0; timeInterval < dr::misc::TimeSteps<Config>; ++timeInterval) {
    timeKernel_.evaluate(
        &coeffs[timeInterval * ConvergenceOrder], timeDerivativePlus, degreesOfFreedomPlus);
    timeKernel_.evaluate(
        &coeffs[timeInterval * ConvergenceOrder], timeDerivativeMinus, degreesOfFreedomMinus);

    const real* plusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps<Config>)
                                   ? &qInterpolatedPlus[timeInterval + 1][0]
                                   : timeDerivativePlusPrefetch;
    const real* minusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps<Config>)
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
    SEISSOL_GPU_PARAM recording::DrConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  // interpolate all timesteps in a single kernel

  runtime.envMany(Cell::NumFaces * dr::misc::NumFaceRelations, [&](void* stream, size_t i) {
    const auto side = i / dr::misc::NumFaceRelations;
    const auto faceRelation = i % dr::misc::NumFaceRelations;

    ConditionalKey minusSideKey(*KernelNames::DrSpaceMap, side, faceRelation);
    if (table.find(minusSideKey) != table.end()) {
      auto& entry = table[minusSideKey];
      const size_t numElements = (entry.get(inner_keys::Dr::Id::IdofsMinus))->getSize();

      auto krnl = gpuCombinedKrnlPrototype_;
      real* tmpMem = reinterpret_cast<real*>(
          device_.api().allocMemAsync(krnl.TmpMaxMemRequiredInBytes * numElements, stream));
      krnl.linearAllocator.initialize(tmpMem);
      krnl.streamPtr = stream;
      krnl.numElements = numElements;

      std::size_t offsetQDR = 0;
      for (std::size_t s = 0; s < dr::misc::TimeSteps<Config>; ++s) {
        krnl.QDR(s) = (entry.get(inner_keys::Dr::Id::QInterpolatedMinus))->getDeviceDataPtr();
        krnl.extraOffset_QDR(s) = offsetQDR;
        offsetQDR += tensor::QDR<Config>::size(s);
      }

      std::size_t offsetDQ = 0;
      for (std::size_t p = 0; p < ConvergenceOrder; ++p) {
        krnl.dQ(p) = const_cast<const real**>(
            (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr());
        krnl.extraOffset_dQ(p) = offsetDQ;
        offsetDQ += tensor::dQ<Config>::size(p);
      }

      for (std::size_t s = 0; s < dr::misc::TimeSteps<Config>; ++s) {
        for (std::size_t p = 0; p < ConvergenceOrder; ++p) {
          krnl.coeffDR(s * ConvergenceOrder + p) = coeffs[s * ConvergenceOrder + p];
        }
      }

      set_I(krnl, entry.get(inner_keys::Dr::Id::IdofsMinus)->getDeviceDataPtr());

      krnl.TinvT =
          const_cast<const real**>((entry.get(inner_keys::Dr::Id::TinvT))->getDeviceDataPtr());
      krnl.execute(side, faceRelation);

      device_.api().freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
    }
  });
#else
  logError() << "No GPU implementation provided";
#endif
}

PerformanceEstimate DynamicRupture::metrics(const DRFaceInformation& faceInfo) const {
  if (isDeviceOn()) {
    return PerformanceEstimate::fromKernel<dynamicRupture::kernel::projectToDR<Config>>(
               faceInfo.plusSide, 0) +
           PerformanceEstimate::fromKernel<dynamicRupture::kernel::projectToDR<Config>>(
               faceInfo.minusSide, faceInfo.faceRelation);
  } else {
    auto estimate = timeKernel_.metrics();

    // 2x evaluateTaylorExpansion
    estimate *= 2;

    estimate += PerformanceEstimate::fromKernel<
        dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Config>>(faceInfo.plusSide,
                                                                                 0);

    estimate += PerformanceEstimate::fromKernel<
        dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Config>>(
        faceInfo.minusSide, faceInfo.faceRelation);

    estimate *= dr::misc::TimeSteps<Config>;

    // legacy CPU memory estimate
    estimate.bytes = (tensor::TinvT<Config>::size() +
                      tensor::QInterpolated<Config>::size() * 2 * dr::misc::TimeSteps<Config> +
                      yateto::computeFamilySize<tensor::dQ<Config>>() * 2) *
                     sizeof(real);

    return estimate;
  }
}

} // namespace seissol::kernels
