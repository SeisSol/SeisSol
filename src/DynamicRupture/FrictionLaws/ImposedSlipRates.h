#ifndef SEISSOL_IMPOSEDSLIPRATES_H
#define SEISSOL_IMPOSEDSLIPRATES_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * Slip rates are set fixed values
 */
template <typename Config, typename STF>
class ImposedSlipRates : public BaseFrictionLaw<Config, ImposedSlipRates<Config, STF>> {
  public:
  using RealT = typename Config::RealT;
  using BaseFrictionLaw<Config, ImposedSlipRates>::BaseFrictionLaw;

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime) {
    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSImposedSlipRates<Config> const* const>(dynRup);
    imposedSlipDirection1 = layerData.var(concreteLts->imposedSlipDirection1);
    imposedSlipDirection2 = layerData.var(concreteLts->imposedSlipDirection2);
    stf.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  void updateFrictionAndSlip(FaultStresses<Config> const& faultStresses,
                             TractionResults<Config>& tractionResults,
                             std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
                             std::array<RealT, misc::numPaddedPoints<Config>>& strengthBuffer,
                             unsigned ltsFace,
                             unsigned timeIndex) {
    const RealT timeIncrement = this->deltaT[timeIndex];
    RealT currentTime = this->mFullUpdateTime;
    for (unsigned i = 0; i <= timeIndex; i++) {
      currentTime += this->deltaT[i];
    }

#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
      const RealT stfEvaluated = stf.evaluate(currentTime, timeIncrement, ltsFace, pointIndex);

      this->traction1[ltsFace][pointIndex] =
          faultStresses.traction1[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * imposedSlipDirection1[ltsFace][pointIndex] * stfEvaluated;
      this->traction2[ltsFace][pointIndex] =
          faultStresses.traction2[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * imposedSlipDirection2[ltsFace][pointIndex] * stfEvaluated;

      this->slipRate1[ltsFace][pointIndex] =
          this->imposedSlipDirection1[ltsFace][pointIndex] * stfEvaluated;
      this->slipRate2[ltsFace][pointIndex] =
          this->imposedSlipDirection2[ltsFace][pointIndex] * stfEvaluated;
      this->slipRateMagnitude[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1[ltsFace][pointIndex], this->slipRate2[ltsFace][pointIndex]);

      // Update slip
      this->slip1[ltsFace][pointIndex] += this->slipRate1[ltsFace][pointIndex] * timeIncrement;
      this->slip2[ltsFace][pointIndex] += this->slipRate2[ltsFace][pointIndex] * timeIncrement;
      this->accumulatedSlipMagnitude[ltsFace][pointIndex] +=
          this->slipRateMagnitude[ltsFace][pointIndex] * timeIncrement;

      tractionResults.traction1[timeIndex][pointIndex] = this->traction1[ltsFace][pointIndex];
      tractionResults.traction2[timeIndex][pointIndex] = this->traction2[ltsFace][pointIndex];
    }
  }

  void preHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
               unsigned ltsFace) {}
  void postHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
                unsigned ltsFace) {}
  void saveDynamicStressOutput(unsigned int ltsFace) {}

  protected:
  RealT (*imposedSlipDirection1)[misc::numPaddedPoints<Config>];
  RealT (*imposedSlipDirection2)[misc::numPaddedPoints<Config>];
  STF stf{};
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_IMPOSEDSLIPRATES_H
