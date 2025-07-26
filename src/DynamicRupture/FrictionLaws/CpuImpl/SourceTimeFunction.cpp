// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SourceTimeFunction.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {
void YoffeSTF::copyLtsTreeToLocal(DynamicRupture::Layer& layerData, real fullUpdateTime) {
  onsetTime = layerData.var<LTSImposedSlipRatesYoffe::OnsetTime>();
  tauS = layerData.var<LTSImposedSlipRatesYoffe::TauS>();
  tauR = layerData.var<LTSImposedSlipRatesYoffe::TauR>();
}

real YoffeSTF::evaluate(real currentTime,
                        [[maybe_unused]] real timeIncrement,
                        size_t ltsFace,
                        uint32_t pointIndex) {
  return regularizedYoffe::regularizedYoffe(currentTime - onsetTime[ltsFace][pointIndex],
                                            tauS[ltsFace][pointIndex],
                                            tauR[ltsFace][pointIndex]);
}

void GaussianSTF::copyLtsTreeToLocal(DynamicRupture::Layer& layerData, real fullUpdateTime) {
  onsetTime = layerData.var<LTSImposedSlipRatesGaussian::OnsetTime>();
  riseTime = layerData.var<LTSImposedSlipRatesGaussian::RiseTime>();
}

real GaussianSTF::evaluate(real currentTime,
                           real timeIncrement,
                           size_t ltsFace,
                           uint32_t pointIndex) {
  const real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
      currentTime - onsetTime[ltsFace][pointIndex], timeIncrement, riseTime[ltsFace][pointIndex]);
  return smoothStepIncrement / timeIncrement;
}

void DeltaSTF::copyLtsTreeToLocal(DynamicRupture::Layer& layerData, real fullUpdateTime) {
  onsetTime = layerData.var<LTSImposedSlipRatesDelta::OnsetTime>();
}

real DeltaSTF::evaluate(real currentTime, real timeIncrement, size_t ltsFace, uint32_t pointIndex) {
  // Currently, the delta pulse is normalized in time equivalent to FL33 and FL34
  return deltaPulse::deltaPulse(currentTime - onsetTime[ltsFace][pointIndex], timeIncrement);
}

} // namespace seissol::dr::friction_law::cpu
