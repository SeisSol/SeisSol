#include "FrictionSolver.h"

namespace seissol::dr::friction_law {

void FrictionSolver::computeDeltaT(double timePoints[CONVERGENCE_ORDER]) {
  deltaT[0] = timePoints[0];
  for (unsigned timeIndex = 1; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }
  // to fill last segment of Gaussian integration
  deltaT[CONVERGENCE_ORDER - 1] = deltaT[CONVERGENCE_ORDER - 1] + deltaT[0];
}

void FrictionSolver::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                        seissol::initializers::DynamicRupture* dynRup,
                                        real fullUpdateTime) {
  impAndEta = layerData.var(dynRup->impAndEta);
  initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS);
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
  mFullUpdateTime = fullUpdateTime;
  averagedSlip = layerData.var(dynRup->averagedSlip);
  dynStressTime = layerData.var(dynRup->dynStressTime);
  dynStressTimePending = layerData.var(dynRup->dynStressTimePending);
  qInterpolatedPlus = layerData.var(dynRup->qInterpolatedPlus);
  qInterpolatedMinus = layerData.var(dynRup->qInterpolatedMinus);
}

real FrictionSolver::calcSmoothStepIncrement(real currentTime, real dt) {
  real gNuc = calcSmoothStep(currentTime);
  real prevTime = currentTime - dt;
  gNuc = gNuc - calcSmoothStep(prevTime);
  return gNuc;
}

real FrictionSolver::calcSmoothStep(real currentTime) {
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < drParameters.t0) {
    return std::exp(misc::power<2>(currentTime - drParameters.t0) /
                    (currentTime * (currentTime - 2.0 * drParameters.t0)));
  } else {
    return 1.0;
  }
}
} // namespace seissol::dr::friction_law
