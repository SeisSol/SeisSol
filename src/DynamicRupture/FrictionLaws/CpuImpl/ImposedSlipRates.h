// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_IMPOSEDSLIPRATES_H_

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law::cpu {
/**
 * Slip rates are set fixed values
 */
template <typename STF>
class ImposedSlipRates : public BaseFrictionLaw<ImposedSlipRates<STF>> {
  public:
  using BaseFrictionLaw<ImposedSlipRates>::BaseFrictionLaw;

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {
    imposedSlipDirection1_ = layerData.var<LTSImposedSlipRates::ImposedSlipDirection1>();
    imposedSlipDirection2_ = layerData.var<LTSImposedSlipRates::ImposedSlipDirection2>();
    stf_.copyStorageToLocal(layerData);
  }

  void updateFrictionAndSlip(const FaultStresses<Executor::Host>& faultStresses,
                             TractionResults<Executor::Host>& tractionResults,
                             std::array<real, misc::NumPaddedPoints>& /*stateVariableBuffer*/,
                             std::array<real, misc::NumPaddedPoints>& /*strengthBuffer*/,
                             std::size_t ltsFace,
                             uint32_t timeIndex) {
    const real timeIncrement = this->deltaT_[timeIndex];
    real currentTime = this->mFullUpdateTime_;
    for (uint32_t i = 0; i <= timeIndex; i++) {
      currentTime += this->deltaT_[i];
    }

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      const real stfEvaluated = stf_.evaluate(currentTime, timeIncrement, ltsFace, pointIndex);

      this->traction1_[ltsFace][pointIndex] = faultStresses.traction1[timeIndex][pointIndex] -
                                              this->impAndEta_[ltsFace].etaS *
                                                  imposedSlipDirection1_[ltsFace][pointIndex] *
                                                  stfEvaluated;
      this->traction2_[ltsFace][pointIndex] = faultStresses.traction2[timeIndex][pointIndex] -
                                              this->impAndEta_[ltsFace].etaS *
                                                  imposedSlipDirection2_[ltsFace][pointIndex] *
                                                  stfEvaluated;

      this->slipRate1_[ltsFace][pointIndex] =
          this->imposedSlipDirection1_[ltsFace][pointIndex] * stfEvaluated;
      this->slipRate2_[ltsFace][pointIndex] =
          this->imposedSlipDirection2_[ltsFace][pointIndex] * stfEvaluated;
      this->slipRateMagnitude_[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1_[ltsFace][pointIndex], this->slipRate2_[ltsFace][pointIndex]);

      // Update slip
      this->slip1_[ltsFace][pointIndex] += this->slipRate1_[ltsFace][pointIndex] * timeIncrement;
      this->slip2_[ltsFace][pointIndex] += this->slipRate2_[ltsFace][pointIndex] * timeIncrement;
      this->accumulatedSlipMagnitude_[ltsFace][pointIndex] +=
          this->slipRateMagnitude_[ltsFace][pointIndex] * timeIncrement;

      tractionResults.traction1[timeIndex][pointIndex] = this->traction1_[ltsFace][pointIndex];
      tractionResults.traction2[timeIndex][pointIndex] = this->traction2_[ltsFace][pointIndex];
    }
  }

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {}
  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {
  }
  void saveDynamicStressOutput(std::size_t ltsFace, real time) {}

  protected:
  real (*__restrict imposedSlipDirection1_)[misc::NumPaddedPoints]{};
  real (*__restrict imposedSlipDirection2_)[misc::NumPaddedPoints]{};
  STF stf_{};
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_IMPOSEDSLIPRATES_H_
