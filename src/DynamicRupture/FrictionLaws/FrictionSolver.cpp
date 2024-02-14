#include "FrictionSolver.h"
#include "generated_code/kernel.h"
#include <yateto/TensorView.h>

namespace seissol::dr::friction_law {

void FrictionSolver::computeDeltaT(const double timePoints[CONVERGENCE_ORDER]) {
  deltaT[0] = timePoints[0];
  for (unsigned timeIndex = 1; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }
  // to fill last segment of Gaussian integration
  deltaT[CONVERGENCE_ORDER - 1] = deltaT[CONVERGENCE_ORDER - 1] + deltaT[0];
}

void FrictionSolver::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                        seissol::initializer::DynamicRupture const* const dynRup,
                                        real fullUpdateTime) {
  impAndEta = layerData.var(dynRup->impAndEta);
  impedanceMatrices = layerData.var(dynRup->impedanceMatrices);
  initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS);
  nucleationStressInFaultCS = layerData.var(dynRup->nucleationStressInFaultCS);
  mu = layerData.var(dynRup->mu);
  accumulatedSlipMagnitude = layerData.var(dynRup->accumulatedSlipMagnitude);
  slip1 = layerData.var(dynRup->slip1);
  slip2 = layerData.var(dynRup->slip2);
  slipRateMagnitude = layerData.var(dynRup->slipRateMagnitude);
  slipRate1 = layerData.var(dynRup->slipRate1);
  slipRate2 = layerData.var(dynRup->slipRate2);
  ruptureTime = layerData.var(dynRup->ruptureTime);
  ruptureTimePending = layerData.var(dynRup->ruptureTimePending);
  peakSlipRate = layerData.var(dynRup->peakSlipRate);
  traction1 = layerData.var(dynRup->traction1);
  traction2 = layerData.var(dynRup->traction2);
  imposedStatePlus = layerData.var(dynRup->imposedStatePlus);
  imposedStateMinus = layerData.var(dynRup->imposedStateMinus);
  energyData = layerData.var(dynRup->drEnergyOutput);
  godunovData = layerData.var(dynRup->godunovData);
  mFullUpdateTime = fullUpdateTime;
  dynStressTime = layerData.var(dynRup->dynStressTime);
  dynStressTimePending = layerData.var(dynRup->dynStressTimePending);
  qInterpolatedPlus = layerData.var(dynRup->qInterpolatedPlus);
  qInterpolatedMinus = layerData.var(dynRup->qInterpolatedMinus);
  initialPressure = layerData.var(dynRup->initialPressure);
  nucleationPressure = layerData.var(dynRup->nucleationPressure);
}
} // namespace seissol::dr::friction_law
