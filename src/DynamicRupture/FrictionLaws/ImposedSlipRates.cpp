#include "ImposedSlipRates.h"

#include "Numerical_aux/RegularizedYoffe.h"

namespace seissol::dr::friction_law {
void ImposedSlipRates::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                          seissol::initializers::DynamicRupture* dynRup,
                                          real fullUpdateTime) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRates*>(dynRup);
  strikeSlip = layerData.var(concreteLts->strikeSlip);
  dipSlip = layerData.var(concreteLts->dipSlip);
  onsetTime = layerData.var(concreteLts->onsetTime);
  tauS = layerData.var(concreteLts->tauS);
  tauR = layerData.var(concreteLts->tauR);
}

void ImposedSlipRates::updateFrictionAndSlip(
    FaultStresses& faultStresses,
    TractionResults& tractionResults,
    std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
    std::array<real, misc::numPaddedPoints>& strengthBuffer,
    unsigned& ltsFace,
    unsigned& timeIndex) {
  real timeIncrement = deltaT[timeIndex];
  real currentTime = mFullUpdateTime;
  for (unsigned i = 0; i <= timeIndex; i++) {
    currentTime += deltaT[i];
  }

  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    // TODO: FL34 with Gauss
    // real gNuc = this->calcSmoothStepIncrement(tn, timeInc) / timeInc;
    real stfEvaluated = regularizedYoffe::regularizedYoffe(currentTime - onsetTime[ltsFace][pointIndex],
                                         tauS[ltsFace][pointIndex],
                                         tauR[ltsFace][pointIndex]);

    tractionResults.traction1[timeIndex][pointIndex] =
        faultStresses.traction1[timeIndex][pointIndex] -
        this->impAndEta[ltsFace].etaS * strikeSlip[ltsFace][pointIndex] * stfEvaluated;
    tractionResults.traction2[timeIndex][pointIndex] =
        faultStresses.traction2[timeIndex][pointIndex] -
        this->impAndEta[ltsFace].etaS * dipSlip[ltsFace][pointIndex] * stfEvaluated;

    this->slipRate1[ltsFace][pointIndex] = this->strikeSlip[ltsFace][pointIndex] * stfEvaluated;
    this->slipRate2[ltsFace][pointIndex] = dipSlip[ltsFace][pointIndex] * stfEvaluated;
    this->slipRateMagnitude[ltsFace][pointIndex] =
        misc::magnitude(this->slipRate1[ltsFace][pointIndex], this->slipRate2[ltsFace][pointIndex]);

    // Update slip
    this->slip1[ltsFace][pointIndex] += this->slipRate1[ltsFace][pointIndex] * timeIncrement;
    this->slip2[ltsFace][pointIndex] += this->slipRate2[ltsFace][pointIndex] * timeIncrement;
    this->accumulatedSlipMagnitude[ltsFace][pointIndex] +=
        this->slipRateMagnitude[ltsFace][pointIndex] * timeIncrement;

    this->traction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
    this->traction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
  }
}

void ImposedSlipRates::preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                               unsigned ltsFace) {}
void ImposedSlipRates::postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                                unsigned ltsFace) {}
void ImposedSlipRates::saveDynamicStressOutput(unsigned int ltsFace) {}

} // namespace seissol::dr::friction_law
