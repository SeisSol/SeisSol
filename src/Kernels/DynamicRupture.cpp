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
#include "Kernels/Common.h"
#include "Kernels/FaceTraces.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
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
#endif

#ifndef NDEBUG
#include <cstdint>
#endif

GENERATE_HAS_MEMBER(faceTraceNeighbor)

namespace seissol::kernels {

// The lift of the fault into a cell is indexed by the side and the face relation. Relation 0
// addresses the plus side, relation 1 the minus side at a zero face orientation index, which the
// canonical vertex numbering guarantees on every interior face.
static_assert(std::size(dynamicRupture::kernel::nodalFlux::ExecutePtrs) ==
              Cell::NumFaces * dr::misc::NumFaceRelations);
static_assert(std::size(tensor::V3mTo2nTWDivM::Size) ==
              Cell::NumFaces * dr::misc::NumFaceRelations);

#ifdef ACL_DEVICE
static_assert(*seissol::recording::DrFaceRelations::Count ==
              Cell::NumFaces * dr::misc::NumFaceRelations);
#endif

void DynamicRupture::setGlobalData(const CompoundGlobalData& global) {
  krnlPrototype_.bindGlobals(*global.onHost);
#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  gpuCombinedKrnlPrototype_.bindGlobals(*global.onDevice);
#endif

  timeKernel_.setGlobalData(global);
}

void DynamicRupture::spaceTimeInterpolation(
    const DRGodunovData* godunovData,
    const real* timeDerivativePlus,
    const real* timeDerivativeMinus,
    real qInterpolatedPlus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
    real qInterpolatedMinus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
    const real* coeffs) {

  assert(timeDerivativePlus != nullptr);
  assert(timeDerivativeMinus != nullptr);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedPlus[0])) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(&qInterpolatedMinus[0])) % Vectorsize == 0);

  // the trace of a side at one point in time
  alignas(PagesizeStack) real tracePlus[tensor::faceTraceNeighbor::size()];
  alignas(PagesizeStack) real traceMinus[tensor::faceTraceNeighbor::size()];

  dynamicRupture::kernel::evaluateAndRotateTraceAtInterpolationPoints krnl = krnlPrototype_;
  kernels::bindFaultRotation(krnl, godunovData->dataTinvT);
  for (std::size_t timeInterval = 0; timeInterval < dr::misc::TimeSteps; ++timeInterval) {
    timeKernel_.evaluateTrace(
        &coeffs[timeInterval * ConvergenceOrder], timeDerivativePlus, tracePlus);
    timeKernel_.evaluateTrace(
        &coeffs[timeInterval * ConvergenceOrder], timeDerivativeMinus, traceMinus);

    krnl.QInterpolated = &qInterpolatedPlus[timeInterval][0];
    krnl.faceTraceNeighbor = tracePlus;
    krnl.execute();

    krnl.QInterpolated = &qInterpolatedMinus[timeInterval][0];
    krnl.faceTraceNeighbor = traceMinus;
    krnl.execute();
  }
}

void DynamicRupture::batchedSpaceTimeInterpolation(
    SEISSOL_GPU_PARAM recording::DrConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  // all points in time of both sides of every face in one kernel: the two sides write their traces
  // in the same parametrisation of the face
  const ConditionalKey key(*KernelNames::DrSpaceMap);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    const auto numElements = (entry.get(inner_keys::Dr::Id::QInterpolated))->getSize();

    auto krnl = gpuCombinedKrnlPrototype_;
    auto tmpMem =
        runtime.memoryHandle<real>((krnl.TmpMaxMemRequiredInBytes * numElements) / sizeof(real));
    krnl.linearAllocator.initialize(tmpMem.get());
    krnl.streamPtr = runtime.stream();
    krnl.numElements = numElements;

    std::size_t offsetQDR = 0;
    for (std::size_t s = 0; s < dr::misc::TimeSteps; ++s) {
      krnl.QDR(s) = (entry.get(inner_keys::Dr::Id::QInterpolated))->getDeviceDataPtr();
      krnl.extraOffset_QDR(s) = offsetQDR;
      offsetQDR += tensor::QDR::size(s);
    }

    kernels::bindFaceDerivativesBatched(
        krnl,
        const_cast<const real**>(
            (entry.get(inner_keys::Dr::Id::FaceDerivatives))->getDeviceDataPtr()));

    for (std::size_t s = 0; s < dr::misc::TimeSteps; ++s) {
      for (std::size_t p = 0; p < ConvergenceOrder; ++p) {
        krnl.coeffDR(s * ConvergenceOrder + p) = coeffs[s * ConvergenceOrder + p];
      }
    }

    // where the kernel forms the trace at a point in time apart
    set_faceTraceNeighbor(krnl, entry.get(inner_keys::Dr::Id::Trace)->getDeviceDataPtr());

    krnl.TinvT =
        const_cast<const real**>((entry.get(inner_keys::Dr::Id::TinvT))->getDeviceDataPtr());
    krnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

PerformanceEstimate DynamicRupture::metrics() const {
  if (isDeviceOn()) {
    return PerformanceEstimate::fromKernel<dynamicRupture::kernel::projectTracesToDR>() * 2;
  } else {
    auto estimate = PerformanceEstimate::fromKernel<kernel::faceDerivativesAtTime>();
    estimate += PerformanceEstimate::fromKernel<
        dynamicRupture::kernel::evaluateAndRotateTraceAtInterpolationPoints>();

    // both sides, every point in time
    estimate *= 2 * dr::misc::TimeSteps;

    // legacy CPU memory estimate
    estimate.bytes = (sizeof(DRGodunovData::dataTinvT) / sizeof(real) +
                      tensor::QInterpolated::size() * 2 * dr::misc::TimeSteps +
                      kernels::Solver::FaceDerivativesSize * 2) *
                     sizeof(real);

    return estimate;
  }
}

} // namespace seissol::kernels
