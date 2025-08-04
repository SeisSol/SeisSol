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
  std::vector<double> deltaT(timePoints.size());
  double sumDt = 0;

  deltaT[0] = timePoints[0];
  sumDt = deltaT[0];
  for (std::size_t timeIndex = 1; timeIndex < timePoints.size(); ++timeIndex) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
    sumDt += deltaT[timeIndex];
  }
  // to fill last segment of Gaussian integration
  deltaT[timePoints.size() - 1] = deltaT[timePoints.size() - 1] + deltaT[0];
  sumDt += deltaT[0];

  return {sumDt, deltaT};
}

template <typename Cfg>
void FrictionSolverImpl<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  const seissol::initializer::AllocationPlace place = allocationPlace();
  impAndEta = layerData.var<DynamicRupture::ImpAndEta>(Cfg(), place);
  impedanceMatrices = layerData.var<DynamicRupture::ImpedanceMatrices>(Cfg(), place);
  initialStressInFaultCS = layerData.var<DynamicRupture::InitialStressInFaultCS>(Cfg(), place);
  nucleationStressInFaultCS =
      layerData.var<DynamicRupture::NucleationStressInFaultCS>(Cfg(), place);
  mu = layerData.var<DynamicRupture::Mu>(Cfg(), place);
  accumulatedSlipMagnitude = layerData.var<DynamicRupture::AccumulatedSlipMagnitude>(Cfg(), place);
  slip1 = layerData.var<DynamicRupture::Slip1>(Cfg(), place);
  slip2 = layerData.var<DynamicRupture::Slip2>(Cfg(), place);
  slipRateMagnitude = layerData.var<DynamicRupture::SlipRateMagnitude>(Cfg(), place);
  slipRate1 = layerData.var<DynamicRupture::SlipRate1>(Cfg(), place);
  slipRate2 = layerData.var<DynamicRupture::SlipRate2>(Cfg(), place);
  ruptureTime = layerData.var<DynamicRupture::RuptureTime>(Cfg(), place);
  ruptureTimePending = layerData.var<DynamicRupture::RuptureTimePending>(Cfg(), place);
  peakSlipRate = layerData.var<DynamicRupture::PeakSlipRate>(Cfg(), place);
  traction1 = layerData.var<DynamicRupture::Traction1>(Cfg(), place);
  traction2 = layerData.var<DynamicRupture::Traction2>(Cfg(), place);
  imposedStatePlus = layerData.var<DynamicRupture::ImposedStatePlus>(Cfg(), place);
  imposedStateMinus = layerData.var<DynamicRupture::ImposedStateMinus>(Cfg(), place);
  energyData = layerData.var<DynamicRupture::DREnergyOutputVar>(Cfg(), place);
  godunovData = layerData.var<DynamicRupture::GodunovData>(Cfg(), place);
  dynStressTime = layerData.var<DynamicRupture::DynStressTime>(Cfg(), place);
  dynStressTimePending = layerData.var<DynamicRupture::DynStressTimePending>(Cfg(), place);
  qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>(Cfg(), place);
  qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>(Cfg(), place);
  initialPressure = layerData.var<DynamicRupture::InitialPressure>(Cfg(), place);
  nucleationPressure = layerData.var<DynamicRupture::NucleationPressure>(Cfg(), place);
}

template <typename Cfg>
seissol::initializer::AllocationPlace FrictionSolverImpl<Cfg>::allocationPlace() {
  return seissol::initializer::AllocationPlace::Host;
}

#define _H_(cfg) template class FrictionSolverImpl<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::dr::friction_law
