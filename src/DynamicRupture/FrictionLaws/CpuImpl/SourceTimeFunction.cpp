// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SourceTimeFunction.h"

#include "Config.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"

#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {
void YoffeSTF::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  onsetTime_ = layerData.var<LTSImposedSlipRatesYoffe::OnsetTime>(Config());
  tauS_ = layerData.var<LTSImposedSlipRatesYoffe::TauS>(Config());
  tauR_ = layerData.var<LTSImposedSlipRatesYoffe::TauR>(Config());
}

real YoffeSTF::evaluate(real currentTime,
                        [[maybe_unused]] real timeIncrement,
                        size_t ltsFace,
                        uint32_t pointIndex) {
  return regularizedYoffe::regularizedYoffe(currentTime - onsetTime_[ltsFace][pointIndex],
                                            tauS_[ltsFace][pointIndex],
                                            tauR_[ltsFace][pointIndex]);
}

void GaussianSTF::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  onsetTime_ = layerData.var<LTSImposedSlipRatesGaussian::OnsetTime>(Config());
  riseTime_ = layerData.var<LTSImposedSlipRatesGaussian::RiseTime>(Config());
}

real GaussianSTF::evaluate(real currentTime,
                           real timeIncrement,
                           size_t ltsFace,
                           uint32_t pointIndex) {
  const real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
      currentTime - onsetTime_[ltsFace][pointIndex], timeIncrement, riseTime_[ltsFace][pointIndex]);
  return smoothStepIncrement / timeIncrement;
}

void DeltaSTF::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  onsetTime_ = layerData.var<LTSImposedSlipRatesDelta::OnsetTime>(Config());
}

real DeltaSTF::evaluate(real currentTime, real timeIncrement, size_t ltsFace, uint32_t pointIndex) {
  // Currently, the delta pulse is normalized in time equivalent to FL33 and FL34
  return deltaPulse::deltaPulse(currentTime - onsetTime_[ltsFace][pointIndex], timeIncrement);
}

} // namespace seissol::dr::friction_law::cpu
