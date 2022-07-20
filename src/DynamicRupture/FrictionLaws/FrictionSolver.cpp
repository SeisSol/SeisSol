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
  mFullUpdateTime = fullUpdateTime;
  averagedSlip = layerData.var(dynRup->averagedSlip);
  dynStressTime = layerData.var(dynRup->dynStressTime);
  dynStressTimePending = layerData.var(dynRup->dynStressTimePending);
  qInterpolatedPlus = layerData.var(dynRup->qInterpolatedPlus);
  qInterpolatedMinus = layerData.var(dynRup->qInterpolatedMinus);
}

void FrictionSolver::adjustInitialStress(size_t ltsFace, size_t timeIndex) {
  if (this->mFullUpdateTime <= this->drParameters.t0) {
    real gNuc = seissol::gaussianNucleationFunction::smoothStepIncrement(
        this->mFullUpdateTime, this->deltaT[timeIndex], this->drParameters.t0);
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      for (unsigned i = 0; i < 6; i++) {
        this->initialStressInFaultCS[ltsFace][pointIndex][i] +=
            nucleationStressInFaultCS[ltsFace][pointIndex][i] * gNuc;
      }
    }
  }
}
} // namespace seissol::dr::friction_law
