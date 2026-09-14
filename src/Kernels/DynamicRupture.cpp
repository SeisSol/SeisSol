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
#include "Equations/Setup.h" // IWYU pragma: keep
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <cassert>
#include <cstring>
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
    real qInterpolatedPlus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
    real qInterpolatedMinus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
    const real* timeDerivativePlusPrefetch,
    const real* timeDerivativeMinusPrefetch,
    const std::vector<TimeCoefficients>& coeffs) {

  // assert alignments
  assert(timeDerivativePlus != nullptr);
  assert(timeDerivativeMinus != nullptr);
  assert((reinterpret_cast<uintptr_t>(timeDerivativePlus)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivativeMinus)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedPlus[0])) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedMinus[0])) % Alignment == 0);
  // What a fault reads of a cell is what the cell transported, evaluated at a
  // point in time -- so these buffers are of that tensor, and the two are the
  // same size wherever a solver's flux is linear. They used to be of the
  // state, with an assertion that the two sizes agree standing in for saying
  // which one was meant.
  alignas(PagesizeStack) real degreesOfFreedomPlus[tensor::I::size()];
  alignas(PagesizeStack) real degreesOfFreedomMinus[tensor::I::size()];

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints krnl = krnlPrototype_;
  for (std::size_t timeInterval = 0; timeInterval < dr::misc::TimeSteps; ++timeInterval) {
    timeKernel_.evaluate(coeffs[timeInterval], timeDerivativePlus, degreesOfFreedomPlus);
    timeKernel_.evaluate(coeffs[timeInterval], timeDerivativeMinus, degreesOfFreedomMinus);

    const real* plusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps)
                                   ? &qInterpolatedPlus[timeInterval + 1][0]
                                   : timeDerivativePlusPrefetch;
    const real* minusPrefetch = (timeInterval + 1 < dr::misc::TimeSteps)
                                    ? &qInterpolatedMinus[timeInterval + 1][0]
                                    : timeDerivativeMinusPrefetch;

    krnl.QInterpolated = &qInterpolatedPlus[timeInterval][0];
    krnl.I = degreesOfFreedomPlus;
    krnl.TinvT = godunovData->dataTinvT;
    krnl._prefetch.QInterpolated = plusPrefetch;
    krnl.execute(faceInfo.plusSide, 0);

    krnl.QInterpolated = &qInterpolatedMinus[timeInterval][0];
    krnl.I = degreesOfFreedomMinus;
    krnl.TinvT = godunovData->dataTinvT;
    krnl._prefetch.QInterpolated = minusPrefetch;
    krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);
  }
}

void DynamicRupture::batchedSpaceTimeInterpolation(
    SEISSOL_GPU_PARAM recording::DrConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM const std::vector<TimeCoefficients>& coeffs,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  // interpolate all timesteps in a single kernel

  runtime.envMany(16, [&](void* stream, size_t i) {
    const auto side = i / 4;
    const auto faceRelation = i % 4;

    ConditionalKey minusSideKey(*KernelNames::DrSpaceMap, side, faceRelation);
    if (table.find(minusSideKey) != table.end()) {
      auto& entry = table[minusSideKey];
      const size_t numElements = (entry.get(inner_keys::Dr::Id::IdofsMinus))->getSize();

      auto krnl = gpuCombinedKrnlPrototype_;
      real* tmpMem = reinterpret_cast<real*>(
          device_.api->allocMemAsync(krnl.TmpMaxMemRequiredInBytes * numElements, stream));
      krnl.linearAllocator.initialize(tmpMem);
      krnl.streamPtr = stream;
      krnl.numElements = numElements;

      std::size_t offsetQDR = 0;
      for (std::size_t s = 0; s < dr::misc::TimeSteps; ++s) {
        krnl.QDR(s) = (entry.get(inner_keys::Dr::Id::QInterpolatedMinus))->getDeviceDataPtr();
        krnl.extraOffset_QDR(s) = offsetQDR;
        offsetQDR += tensor::QDR::size(s);
      }

      std::size_t offsetDQ = 0;
      for (std::size_t p = 0; p < ConvergenceOrder; ++p) {
        krnl.dQ(p) = const_cast<const real**>(
            (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr());
        krnl.extraOffset_dQ(p) = offsetDQ;
        offsetDQ += tensor::dQ::size(p);
      }

      // What the cell carries beyond its state, behind the expansion of the
      // state in the same buffer. A fault reads both -- the state's columns
      // from the Taylor sum and the rest from this -- and the coefficients for
      // the second are set below whether or not there is anything to apply
      // them to, so leaving this unbound sums one expansion and dereferences
      // nothing for the other.
      std::size_t offsetTransport = yateto::computeFamilySize<tensor::dQ>();
      for (std::size_t p = 0; p < yateto::numFamilyMembers<tensor::transportDer>(); ++p) {
        krnl.transportDer(p) = const_cast<const real**>(
            (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr());
        krnl.extraOffset_transportDer(p) = offsetTransport;
        offsetTransport += tensor::transportDer::size(p);
      }

      for (std::size_t s = 0; s < dr::misc::TimeSteps; ++s) {
        for (std::size_t p = 0; p < ConvergenceOrder; ++p) {
          seissol::model::bindFaultTimeCoefficient<model::MaterialT>(
              krnl, s * ConvergenceOrder + p, coeffs[s], p);
        }
      }

      set_I(krnl, entry.get(inner_keys::Dr::Id::IdofsMinus)->getDeviceDataPtr());

      krnl.TinvT =
          const_cast<const real**>((entry.get(inner_keys::Dr::Id::TinvT))->getDeviceDataPtr());
      krnl.execute(side, faceRelation);

      device_.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
    }
  });
#else
  logError() << "No GPU implementation provided";
#endif
}

PerformanceEstimate DynamicRupture::metrics(const DRFaceInformation& faceInfo) const {
  if (isDeviceOn()) {
    return PerformanceEstimate::fromKernel<dynamicRupture::kernel::projectToDR>(faceInfo.plusSide,
                                                                                0) +
           PerformanceEstimate::fromKernel<dynamicRupture::kernel::projectToDR>(
               faceInfo.minusSide, faceInfo.faceRelation);
  } else {
    auto estimate = timeKernel_.metrics();

    // 2x evaluateTaylorExpansion
    estimate *= 2;

    estimate += PerformanceEstimate::fromKernel<
        dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints>(faceInfo.plusSide, 0);

    estimate += PerformanceEstimate::fromKernel<
        dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints>(faceInfo.minusSide,
                                                                         faceInfo.faceRelation);

    estimate *= dr::misc::TimeSteps;

    // legacy CPU memory estimate. Both sides of the face, and both expansions
    // per side where a cell keeps two: a fault reads the state's columns from
    // the one and everything else from the other.
    estimate.bytes =
        (tensor::TinvT::size() + tensor::QInterpolated::size() * 2 * dr::misc::TimeSteps +
         (yateto::computeFamilySize<tensor::dQ>() +
          yateto::computeFamilySize<tensor::transportDer>()) *
             2) *
        sizeof(real);

    return estimate;
  }
}

} // namespace seissol::kernels
