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

void FrictionSolver::copyLtsTreeToLocal(DynamicRupture::Layer& layerData, real fullUpdateTime) {
  const seissol::initializer::AllocationPlace place = allocationPlace();
  impAndEta = layerData.var<DynamicRupture::ImpAndEta>(place);
  impedanceMatrices = layerData.var<DynamicRupture::ImpedanceMatrices>(place);
  initialStressInFaultCS = layerData.var<DynamicRupture::InitialStressInFaultCS>(place);
  nucleationStressInFaultCS = layerData.var<DynamicRupture::NucleationStressInFaultCS>(place);
  mu = layerData.var<DynamicRupture::Mu>(place);
  accumulatedSlipMagnitude = layerData.var<DynamicRupture::AccumulatedSlipMagnitude>(place);
  slip1 = layerData.var<DynamicRupture::Slip1>(place);
  slip2 = layerData.var<DynamicRupture::Slip2>(place);
  slipRateMagnitude = layerData.var<DynamicRupture::SlipRateMagnitude>(place);
  slipRate1 = layerData.var<DynamicRupture::SlipRate1>(place);
  slipRate2 = layerData.var<DynamicRupture::SlipRate2>(place);
  ruptureTime = layerData.var<DynamicRupture::RuptureTime>(place);
  ruptureTimePending = layerData.var<DynamicRupture::RuptureTimePending>(place);
  peakSlipRate = layerData.var<DynamicRupture::PeakSlipRate>(place);
  traction1 = layerData.var<DynamicRupture::Traction1>(place);
  traction2 = layerData.var<DynamicRupture::Traction2>(place);
  imposedStatePlus = layerData.var<DynamicRupture::ImposedStatePlus>(place);
  imposedStateMinus = layerData.var<DynamicRupture::ImposedStateMinus>(place);
  energyData = layerData.var<DynamicRupture::DREnergyOutputVar>(place);
  godunovData = layerData.var<DynamicRupture::GodunovData>(place);
  mFullUpdateTime = fullUpdateTime;
  dynStressTime = layerData.var<DynamicRupture::DynStressTime>(place);
  dynStressTimePending = layerData.var<DynamicRupture::DynStressTimePending>(place);
  qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>(place);
  qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>(place);
  initialPressure = layerData.var<DynamicRupture::InitialPressure>(place);
  nucleationPressure = layerData.var<DynamicRupture::NucleationPressure>(place);
}
seissol::initializer::AllocationPlace FrictionSolver::allocationPlace() {
  return seissol::initializer::AllocationPlace::Host;
}

} // namespace seissol::dr::friction_law
