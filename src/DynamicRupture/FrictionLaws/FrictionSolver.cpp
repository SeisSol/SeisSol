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

template <typename Cfg>
FrictionSolver::FrictionTime FrictionSolver::computeDeltaT(const std::vector<double>& timePoints) {
  std::vector<double> deltaT(misc::TimeSteps<Cfg>);

  if (timePoints.size() != deltaT.size()) {
    logError() << "Internal time point count mismatch. Given vs. expected:" << timePoints.size()
               << deltaT.size();
  }

  deltaT[0] = timePoints[0]; // - 0
  for (std::size_t timeIndex = 1; timeIndex < misc::TimeSteps<Cfg>; ++timeIndex) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }

  // add the segment [lastPoint, timestep] to the last point
  deltaT.back() += timePoints[0];

  return {deltaT};
}

template <typename Cfg>
void FrictionSolverImpl<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  const seissol::initializer::AllocationPlace place = allocationPlace();
  impAndEta_ = layerData.var<DynamicRupture::ImpAndEta>(Cfg(), place);
  impedanceMatrices_ = layerData.var<DynamicRupture::ImpedanceMatrices>(Cfg(), place);
  stressSourceInFaultCS_ = layerData.var<DynamicRupture::StressSourceInFaultCS>(Cfg(), place);
  mu_ = layerData.var<DynamicRupture::Mu>(Cfg(), place);
  accumulatedSlipMagnitude_ = layerData.var<DynamicRupture::AccumulatedSlipMagnitude>(Cfg(), place);
  slip1_ = layerData.var<DynamicRupture::Slip1>(Cfg(), place);
  slip2_ = layerData.var<DynamicRupture::Slip2>(Cfg(), place);
  slipRateMagnitude_ = layerData.var<DynamicRupture::SlipRateMagnitude>(Cfg(), place);
  slipRate1_ = layerData.var<DynamicRupture::SlipRate1>(Cfg(), place);
  slipRate2_ = layerData.var<DynamicRupture::SlipRate2>(Cfg(), place);
  ruptureTime_ = layerData.var<DynamicRupture::RuptureTime>(Cfg(), place);
  ruptureTimePending_ = layerData.var<DynamicRupture::RuptureTimePending>(Cfg(), place);
  peakSlipRate_ = layerData.var<DynamicRupture::PeakSlipRate>(Cfg(), place);
  traction1_ = layerData.var<DynamicRupture::Traction1>(Cfg(), place);
  traction2_ = layerData.var<DynamicRupture::Traction2>(Cfg(), place);
  imposedStatePlus_ = layerData.var<DynamicRupture::ImposedStatePlus>(Cfg(), place);
  imposedStateMinus_ = layerData.var<DynamicRupture::ImposedStateMinus>(Cfg(), place);
  energyData_ = layerData.var<DynamicRupture::DREnergyOutputVar>(Cfg(), place);
  godunovData_ = layerData.var<DynamicRupture::GodunovData>(Cfg(), place);
  dynStressTime_ = layerData.var<DynamicRupture::DynStressTime>(Cfg(), place);
  dynStressTimePending_ = layerData.var<DynamicRupture::DynStressTimePending>(Cfg(), place);
  qInterpolatedPlus_ = layerData.var<DynamicRupture::QInterpolatedPlus>(Cfg(), place);
  qInterpolatedMinus_ = layerData.var<DynamicRupture::QInterpolatedMinus>(Cfg(), place);
  stressSourcePressure_ = layerData.var<DynamicRupture::StressSourcePressure>(Cfg(), place);
  stressSourceOnset_ = layerData.var<DynamicRupture::StressSourceOnset>(Cfg(), place);
  stressSourceRiseTime_ = layerData.var<DynamicRupture::StressSourceRiseTime>(Cfg(), place);
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template FrictionSolver::FrictionTime FrictionSolver::computeDeltaT<Cfg>(                        \
      const std::vector<double>& timePoints);                                                      \
  template class FrictionSolverImpl<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::dr::friction_law
