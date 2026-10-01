// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FrictionSolver.h"

#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"

#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::friction_law {

FrictionSolver::FrictionTime FrictionSolver::computeDeltaT(const std::vector<double>& timePoints) {
  std::vector<double> deltaT(misc::TimeSteps<Config>);

  if (timePoints.size() != deltaT.size()) {
    logError() << "Internal time point count mismatch. Given vs. expected:" << timePoints.size()
               << deltaT.size();
  }

  deltaT[0] = timePoints[0]; // - 0
  for (std::size_t timeIndex = 1; timeIndex < misc::TimeSteps<Config>; ++timeIndex) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }

  // add the segment [lastPoint, timestep] to the last point
  deltaT.back() += timePoints[0];

  return {deltaT};
}

void FrictionSolver::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  const seissol::initializer::AllocationPlace place = allocationPlace();
  impAndEta_ = layerData.var<DynamicRupture::ImpAndEta>(Config(), place);
  impedanceMatrices_ = layerData.var<DynamicRupture::ImpedanceMatrices>(Config(), place);
  stressSourceInFaultCS_ = layerData.var<DynamicRupture::StressSourceInFaultCS>(Config(), place);
  mu_ = layerData.var<DynamicRupture::Mu>(Config(), place);
  accumulatedSlipMagnitude_ =
      layerData.var<DynamicRupture::AccumulatedSlipMagnitude>(Config(), place);
  slip1_ = layerData.var<DynamicRupture::Slip1>(Config(), place);
  slip2_ = layerData.var<DynamicRupture::Slip2>(Config(), place);
  slipRateMagnitude_ = layerData.var<DynamicRupture::SlipRateMagnitude>(Config(), place);
  slipRate1_ = layerData.var<DynamicRupture::SlipRate1>(Config(), place);
  slipRate2_ = layerData.var<DynamicRupture::SlipRate2>(Config(), place);
  ruptureTime_ = layerData.var<DynamicRupture::RuptureTime>(Config(), place);
  ruptureTimePending_ = layerData.var<DynamicRupture::RuptureTimePending>(Config(), place);
  peakSlipRate_ = layerData.var<DynamicRupture::PeakSlipRate>(Config(), place);
  traction1_ = layerData.var<DynamicRupture::Traction1>(Config(), place);
  traction2_ = layerData.var<DynamicRupture::Traction2>(Config(), place);
  imposedStatePlus_ = layerData.var<DynamicRupture::ImposedStatePlus>(Config(), place);
  imposedStateMinus_ = layerData.var<DynamicRupture::ImposedStateMinus>(Config(), place);
  energyData_ = layerData.var<DynamicRupture::DREnergyOutputVar>(Config(), place);
  godunovData_ = layerData.var<DynamicRupture::GodunovData>(Config(), place);
  dynStressTime_ = layerData.var<DynamicRupture::DynStressTime>(Config(), place);
  dynStressTimePending_ = layerData.var<DynamicRupture::DynStressTimePending>(Config(), place);
  qInterpolatedPlus_ = layerData.var<DynamicRupture::QInterpolatedPlus>(Config(), place);
  qInterpolatedMinus_ = layerData.var<DynamicRupture::QInterpolatedMinus>(Config(), place);
  stressSourcePressure_ = layerData.var<DynamicRupture::StressSourcePressure>(Config(), place);
  stressSourceOnset_ = layerData.var<DynamicRupture::StressSourceOnset>(Config(), place);
  stressSourceRiseTime_ = layerData.var<DynamicRupture::StressSourceRiseTime>(Config(), place);
}
} // namespace seissol::dr::friction_law
