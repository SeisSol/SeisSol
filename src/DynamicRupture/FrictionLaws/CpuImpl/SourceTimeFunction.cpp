// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SourceTimeFunction.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"
#include <cstddef>

namespace seissol::dr::friction_law::cpu {
void YoffeSTF::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                  const seissol::initializer::DynamicRupture* dynRup,
                                  real fullUpdateTime) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesYoffe*>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
  tauS = layerData.var(concreteLts->tauS);
  tauR = layerData.var(concreteLts->tauR);
}

real YoffeSTF::evaluate(real currentTime,
                        [[maybe_unused]] real timeIncrement,
                        size_t ltsFace,
                        size_t pointIndex) {
  return regularizedYoffe::regularizedYoffe(currentTime - onsetTime[ltsFace][pointIndex],
                                            tauS[ltsFace][pointIndex],
                                            tauR[ltsFace][pointIndex]);
}

void GaussianSTF::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* dynRup,
                                     real fullUpdateTime) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesGaussian*>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
  riseTime = layerData.var(concreteLts->riseTime);
}

real GaussianSTF::evaluate(real currentTime,
                           real timeIncrement,
                           size_t ltsFace,
                           size_t pointIndex) {
  const real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
      currentTime - onsetTime[ltsFace][pointIndex], timeIncrement, riseTime[ltsFace][pointIndex]);
  return smoothStepIncrement / timeIncrement;
}

void DeltaSTF::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                  const seissol::initializer::DynamicRupture* const dynRup,
                                  real fullUpdateTime) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesDelta*>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
}

real DeltaSTF::evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex) {
  // Currently, the delta pulse is normalized in time equivalent to FL33 and FL34
  return deltaPulse::deltaPulse(currentTime - onsetTime[ltsFace][pointIndex], timeIncrement);
}

} // namespace seissol::dr::friction_law::cpu
