// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FrictionSolver.h"

#include "Common/Constants.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include <cstddef>
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

  // use that time points are symmetric to compute dt
  const auto sumDt = timePoints.back() + timePoints[0];

  return {sumDt, deltaT};
}

void FrictionSolver::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                        const seissol::initializer::DynamicRupture* const dynRup,
                                        real fullUpdateTime) {
  const seissol::initializer::AllocationPlace place = allocationPlace();
  impAndEta = layerData.var(dynRup->impAndEta, place);
  impedanceMatrices = layerData.var(dynRup->impedanceMatrices, place);
  initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS, place);
  for (std::size_t i = 0; i < dynRup->nucleationStressInFaultCS.size(); ++i) {
    nucleationStressInFaultCS[i] = layerData.var(dynRup->nucleationStressInFaultCS[i], place);
  }
  mu = layerData.var(dynRup->mu, place);
  accumulatedSlipMagnitude = layerData.var(dynRup->accumulatedSlipMagnitude, place);
  slip1 = layerData.var(dynRup->slip1, place);
  slip2 = layerData.var(dynRup->slip2, place);
  slipRateMagnitude = layerData.var(dynRup->slipRateMagnitude, place);
  slipRate1 = layerData.var(dynRup->slipRate1, place);
  slipRate2 = layerData.var(dynRup->slipRate2, place);
  ruptureTime = layerData.var(dynRup->ruptureTime, place);
  ruptureTimePending = layerData.var(dynRup->ruptureTimePending, place);
  peakSlipRate = layerData.var(dynRup->peakSlipRate, place);
  traction1 = layerData.var(dynRup->traction1, place);
  traction2 = layerData.var(dynRup->traction2, place);
  imposedStatePlus = layerData.var(dynRup->imposedStatePlus, place);
  imposedStateMinus = layerData.var(dynRup->imposedStateMinus, place);
  energyData = layerData.var(dynRup->drEnergyOutput, place);
  godunovData = layerData.var(dynRup->godunovData, place);
  mFullUpdateTime = fullUpdateTime;
  dynStressTime = layerData.var(dynRup->dynStressTime, place);
  dynStressTimePending = layerData.var(dynRup->dynStressTimePending, place);
  qInterpolatedPlus = layerData.var(dynRup->qInterpolatedPlus, place);
  qInterpolatedMinus = layerData.var(dynRup->qInterpolatedMinus, place);
  initialPressure = layerData.var(dynRup->initialPressure, place);
  for (std::size_t i = 0; i < dynRup->nucleationPressure.size(); ++i) {
    nucleationPressure[i] = layerData.var(dynRup->nucleationPressure[i], place);
  }
}
seissol::initializer::AllocationPlace FrictionSolver::allocationPlace() {
  return seissol::initializer::AllocationPlace::Host;
}

} // namespace seissol::dr::friction_law
