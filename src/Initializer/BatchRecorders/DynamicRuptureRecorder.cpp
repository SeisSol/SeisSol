// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"

#include <cstddef>
#include <vector>

using namespace seissol::initializer;
using namespace seissol::recording;

void DynamicRuptureRecorder::record(DynamicRupture::Layer& layer) {
  setUpContext(layer);
  recordSpaceInterpolation();
}

void DynamicRuptureRecorder::recordSpaceInterpolation() {
  auto* qInterpolatedPlus =
      currentLayer_->var<DynamicRupture::QInterpolatedPlus>(AllocationPlace::Device);
  auto* qInterpolatedMinus =
      currentLayer_->var<DynamicRupture::QInterpolatedMinus>(AllocationPlace::Device);
  real* tracesScratch = static_cast<real*>(
      currentLayer_->var<DynamicRupture::TracesOnDevice>(AllocationPlace::Device));

  DRGodunovData* godunovData =
      currentLayer_->var<DynamicRupture::GodunovData>(AllocationPlace::Device);

  real* const* timeDerivativePlus = currentLayer_->var<DynamicRupture::TimeDerivativePlusDevice>();
  real* const* timeDerivativeMinus =
      currentLayer_->var<DynamicRupture::TimeDerivativeMinusDevice>();

  const auto size = currentLayer_->size();
  if (size > 0) {
    // Both sides of every face at once: they write their traces in the same parametrisation of
    // the face, which one kernel reads for either.
    std::vector<real*> faceDerivativesPtrs;
    std::vector<real*> tracePtrs;
    std::vector<real*> qInterpolatedPtrs;
    std::vector<real*> tInvTPtrs;
    faceDerivativesPtrs.reserve(2 * size);
    tracePtrs.reserve(2 * size);
    qInterpolatedPtrs.reserve(2 * size);
    tInvTPtrs.reserve(2 * size);

    for (std::size_t faceId = 0; faceId < size; ++faceId) {
      faceDerivativesPtrs.push_back(timeDerivativePlus[faceId]);
      tracePtrs.push_back(tracesScratch + 2 * faceId * kernels::Solver::FaceTraceSize);
      qInterpolatedPtrs.push_back(&qInterpolatedPlus[faceId][0][0]);
      tInvTPtrs.push_back(godunovData[faceId].dataTinvT);

      faceDerivativesPtrs.push_back(timeDerivativeMinus[faceId]);
      tracePtrs.push_back(tracesScratch + (2 * faceId + 1) * kernels::Solver::FaceTraceSize);
      qInterpolatedPtrs.push_back(&qInterpolatedMinus[faceId][0][0]);
      tInvTPtrs.push_back(godunovData[faceId].dataTinvT);
    }

    const ConditionalKey key(*KernelNames::DrSpaceMap);
    (*currentDrTable_)[key].set(inner_keys::Dr::Id::FaceDerivatives, faceDerivativesPtrs);
    (*currentDrTable_)[key].set(inner_keys::Dr::Id::Trace, tracePtrs);
    (*currentDrTable_)[key].set(inner_keys::Dr::Id::QInterpolated, qInterpolatedPtrs);
    (*currentDrTable_)[key].set(inner_keys::Dr::Id::TinvT, tInvTPtrs);
  }
}
