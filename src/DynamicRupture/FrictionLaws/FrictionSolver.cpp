// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FrictionSolver.h"

#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"

#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::friction_law {

FrictionSolver::FrictionTime FrictionSolver::computeDeltaT(const std::vector<double>& timePoints) {
  std::vector<double> deltaT(misc::TimeSteps);

  if (timePoints.size() != deltaT.size()) {
    logError() << "Internal time point count mismatch. Given vs. expected:" << timePoints.size()
               << deltaT.size();
  }

  deltaT[0] = timePoints[0]; // - 0
  for (std::size_t timeIndex = 1; timeIndex < misc::TimeSteps; ++timeIndex) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }

  // add the segment [lastPoint, timestep] to the last point
  deltaT.back() += timePoints[0];

  return {deltaT};
}

void FrictionSolver::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  const seissol::initializer::AllocationPlace place = allocationPlace();
  impAndEta_ = layerData.var<DynamicRupture::ImpAndEta>(place);
  impedanceMatrices_ = layerData.var<DynamicRupture::ImpedanceMatrices>(place);
  initialStressInFaultCS_ = layerData.var<DynamicRupture::InitialStressInFaultCS>(place);
  nucleationStressInFaultCS_ = layerData.var<DynamicRupture::NucleationStressInFaultCS>(place);
  mu_ = layerData.var<DynamicRupture::Mu>(place);
  accumulatedSlipMagnitude_ = layerData.var<DynamicRupture::AccumulatedSlipMagnitude>(place);
  slip1_ = layerData.var<DynamicRupture::Slip1>(place);
  slip2_ = layerData.var<DynamicRupture::Slip2>(place);
  slipRateMagnitude_ = layerData.var<DynamicRupture::SlipRateMagnitude>(place);
  slipRate1_ = layerData.var<DynamicRupture::SlipRate1>(place);
  slipRate2_ = layerData.var<DynamicRupture::SlipRate2>(place);
  ruptureTime_ = layerData.var<DynamicRupture::RuptureTime>(place);
  ruptureTimePending_ = layerData.var<DynamicRupture::RuptureTimePending>(place);
  peakSlipRate_ = layerData.var<DynamicRupture::PeakSlipRate>(place);
  traction1_ = layerData.var<DynamicRupture::Traction1>(place);
  traction2_ = layerData.var<DynamicRupture::Traction2>(place);
  imposedStatePlus_ = layerData.var<DynamicRupture::ImposedStatePlus>(place);
  imposedStateMinus_ = layerData.var<DynamicRupture::ImposedStateMinus>(place);
  energyData_ = layerData.var<DynamicRupture::DREnergyOutputVar>(place);
  godunovData_ = layerData.var<DynamicRupture::GodunovData>(place);
  dynStressTime_ = layerData.var<DynamicRupture::DynStressTime>(place);
  dynStressTimePending_ = layerData.var<DynamicRupture::DynStressTimePending>(place);
  qInterpolatedPlus_ = layerData.var<DynamicRupture::QInterpolatedPlus>(place);
  qInterpolatedMinus_ = layerData.var<DynamicRupture::QInterpolatedMinus>(place);
  initialPressure_ = layerData.var<DynamicRupture::InitialPressure>(place);
  nucleationPressure_ = layerData.var<DynamicRupture::NucleationPressure>(place);
}
seissol::initializer::AllocationPlace FrictionSolver::allocationPlace() {
  return seissol::initializer::AllocationPlace::Host;
}

} // namespace seissol::dr::friction_law
