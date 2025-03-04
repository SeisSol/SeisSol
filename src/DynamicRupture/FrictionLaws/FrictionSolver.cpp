// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
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

namespace seissol::dr::friction_law {

void FrictionSolver::computeDeltaT(const double timePoints[ConvergenceOrder]) {
  deltaT[0] = timePoints[0];
  sumDt = deltaT[0];
  for (std::size_t timeIndex = 1; timeIndex < ConvergenceOrder; timeIndex++) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
    sumDt += deltaT[timeIndex];
  }
  // to fill last segment of Gaussian integration
  deltaT[ConvergenceOrder - 1] = deltaT[ConvergenceOrder - 1] + deltaT[0];
  sumDt += deltaT[0];
}

void FrictionSolver::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                        const seissol::initializer::DynamicRupture* const dynRup,
                                        real fullUpdateTime) {
  const seissol::initializer::AllocationPlace place = allocationPlace();
  impAndEta = layerData.var(dynRup->impAndEta, place);
  impedanceMatrices = layerData.var(dynRup->impedanceMatrices, place);
  initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS, place);
  nucleationStressInFaultCS = layerData.var(dynRup->nucleationStressInFaultCS, place);
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
  nucleationPressure = layerData.var(dynRup->nucleationPressure, place);
}
seissol::initializer::AllocationPlace FrictionSolver::allocationPlace() {
  return seissol::initializer::AllocationPlace::Host;
}

} // namespace seissol::dr::friction_law
