#ifndef SEISSOL_IMPOSEDSLIPRATES_H
#define SEISSOL_IMPOSEDSLIPRATES_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * Slip rates are set fixed values
 */
template <typename STF>
class ImposedSlipRates : public BaseFrictionLaw<ImposedSlipRates<STF>> {
  public:
  using BaseFrictionLaw<ImposedSlipRates>::BaseFrictionLaw;

  real (*strikeSlip)[misc::numPaddedPoints];
  real (*dipSlip)[misc::numPaddedPoints];

  STF stf{};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRates*>(dynRup);
    strikeSlip = layerData.var(concreteLts->strikeSlip);
    dipSlip = layerData.var(concreteLts->dipSlip);
    stf.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             TractionResults& tractionResults,
                             std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::numPaddedPoints>& strengthBuffer,
                             unsigned& ltsFace,
                             unsigned& timeIndex) {
    real timeIncrement = this->deltaT[timeIndex];
    real currentTime = this->mFullUpdateTime;
    for (unsigned i = 0; i <= timeIndex; i++) {
      currentTime += this->deltaT[i];
    }

    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      real stfEvaluated = stf.evaluate(currentTime, timeIncrement, ltsFace, pointIndex);

      tractionResults.traction1[timeIndex][pointIndex] =
          faultStresses.xyStress[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * strikeSlip[ltsFace][pointIndex] * stfEvaluated;
      tractionResults.traction2[timeIndex][pointIndex] =
          faultStresses.xzStress[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * dipSlip[ltsFace][pointIndex] * stfEvaluated;

      this->slipRate1[ltsFace][pointIndex] = this->strikeSlip[ltsFace][pointIndex] * stfEvaluated;
      this->slipRate2[ltsFace][pointIndex] = this->dipSlip[ltsFace][pointIndex] * stfEvaluated;
      this->slipRateMagnitude[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1[ltsFace][pointIndex], this->slipRate2[ltsFace][pointIndex]);

      // Update slip
      this->slip1[ltsFace][pointIndex] += this->slipRate1[ltsFace][pointIndex] * timeIncrement;
      this->slip2[ltsFace][pointIndex] += this->slipRate2[ltsFace][pointIndex] * timeIncrement;
      this->accumulatedSlipMagnitude[ltsFace][pointIndex] +=
          this->slipRateMagnitude[ltsFace][pointIndex] * timeIncrement;

      this->tractionXY[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
      this->tractionXZ[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
    }
  }

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {}
  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {}
  void saveDynamicStressOutput(unsigned int ltsFace) {}
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_IMPOSEDSLIPRATES_H
