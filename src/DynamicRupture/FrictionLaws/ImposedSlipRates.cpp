#include "ImposedSlipRates.h"

namespace seissol::dr::friction_law {
void ImposedSlipRates::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                          seissol::initializers::DynamicRupture* dynRup,
                                          real fullUpdateTime) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRates*>(dynRup);
  nucleationStressInFaultCS = layerData.var(concreteLts->nucleationStressInFaultCS);
  averagedSlip = layerData.var(concreteLts->averagedSlip);
}

void ImposedSlipRates::updateFrictionAndSlip(
    FaultStresses& faultStresses,
    TractionResults& tractionResults,
    std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
    std::array<real, misc::numPaddedPoints>& strengthBuffer,
    unsigned& ltsFace,
    unsigned& timeIndex) {
  real timeInc = deltaT[timeIndex];
  real tn = mFullUpdateTime;
  for (unsigned i = 0; i <= timeIndex; i++) {
    tn += deltaT[i];
  }
  real gNuc = this->calcSmoothStepIncrement(tn, timeInc) / timeInc;

  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    //! EQN%NucleationStressInFaultCS (1 and 2) contains the slip in FaultCS
    tractionResults.xyTraction[timeIndex][pointIndex] =
        faultStresses.xyStress[timeIndex][pointIndex] -
        this->impAndEta[ltsFace].etaS * this->nucleationStressInFaultCS[ltsFace][pointIndex][0] *
            gNuc;
    tractionResults.xzTraction[timeIndex][pointIndex] =
        faultStresses.xzStress[timeIndex][pointIndex] -
        this->impAndEta[ltsFace].etaS * this->nucleationStressInFaultCS[ltsFace][pointIndex][1] *
            gNuc;
    this->slipRate1[ltsFace][pointIndex] =
        this->nucleationStressInFaultCS[ltsFace][pointIndex][0] * gNuc;
    this->slipRate2[ltsFace][pointIndex] =
        this->nucleationStressInFaultCS[ltsFace][pointIndex][1] * gNuc;
    this->slipRateMagnitude[ltsFace][pointIndex] =
        misc::magnitude(this->slipRate1[ltsFace][pointIndex], this->slipRate2[ltsFace][pointIndex]);

    //! Update slip
    this->slip1[ltsFace][pointIndex] += this->slipRate1[ltsFace][pointIndex] * timeInc;
    this->slip2[ltsFace][pointIndex] += this->slipRate2[ltsFace][pointIndex] * timeInc;
    this->accumulatedSlipMagnitude[ltsFace][pointIndex] +=
        this->slipRateMagnitude[ltsFace][pointIndex] * timeInc;

    this->tractionXY[ltsFace][pointIndex] = tractionResults.xyTraction[timeIndex][pointIndex];
    this->tractionXZ[ltsFace][pointIndex] = tractionResults.xzTraction[timeIndex][pointIndex];
  }
}

void ImposedSlipRates::preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                               unsigned ltsFace) {}
void ImposedSlipRates::postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                                unsigned ltsFace) {}
void ImposedSlipRates::saveDynamicStressOutput(unsigned int ltsFace) {}
} // namespace seissol::dr::friction_law
