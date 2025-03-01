// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_IMPOSEDSLIPRATES_H_

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law::cpu {
/**
 * Slip rates are set fixed values
 */
template <typename STF>
class ImposedSlipRates : public BaseFrictionLaw<ImposedSlipRates<STF>> {
  public:
  using BaseFrictionLaw<ImposedSlipRates>::BaseFrictionLaw;

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSImposedSlipRates*>(dynRup);
    imposedSlipDirection1 = layerData.var(concreteLts->imposedSlipDirection1);
    imposedSlipDirection2 = layerData.var(concreteLts->imposedSlipDirection2);
    stf.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  void updateFrictionAndSlip(const FaultStresses<Executor::Host>& faultStresses,
                             TractionResults<Executor::Host>& tractionResults,
                             std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::NumPaddedPoints>& strengthBuffer,
                             unsigned ltsFace,
                             unsigned timeIndex) {
    const real timeIncrement = this->deltaT[timeIndex];
    real currentTime = this->mFullUpdateTime;
    for (unsigned i = 0; i <= timeIndex; i++) {
      currentTime += this->deltaT[i];
    }

#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      const real stfEvaluated = stf.evaluate(currentTime, timeIncrement, ltsFace, pointIndex);

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

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {}
  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {}
  void saveDynamicStressOutput(unsigned int ltsFace) {}

  protected:
  real (*imposedSlipDirection1)[misc::NumPaddedPoints]{};
  real (*imposedSlipDirection2)[misc::NumPaddedPoints]{};
  STF stf{};
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_IMPOSEDSLIPRATES_H_
