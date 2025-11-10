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
template <typename Cfg>
void YoffeSTF<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  onsetTime = layerData.var<LTSImposedSlipRatesYoffe::OnsetTime>(Cfg());
  tauS = layerData.var<LTSImposedSlipRatesYoffe::TauS>(Cfg());
  tauR = layerData.var<LTSImposedSlipRatesYoffe::TauR>(Cfg());
}

template <typename Cfg>
Real<Cfg> YoffeSTF<Cfg>::evaluate(real currentTime,
                                  [[maybe_unused]] real timeIncrement,
                                  size_t ltsFace,
                                  uint32_t pointIndex) {
  return regularizedYoffe::regularizedYoffe(currentTime - onsetTime[ltsFace][pointIndex],
                                            tauS[ltsFace][pointIndex],
                                            tauR[ltsFace][pointIndex]);
}

template <typename Cfg>
void GaussianSTF<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  onsetTime = layerData.var<LTSImposedSlipRatesGaussian::OnsetTime>(Cfg());
  riseTime = layerData.var<LTSImposedSlipRatesGaussian::RiseTime>(Cfg());
}

template <typename Cfg>
Real<Cfg> GaussianSTF<Cfg>::evaluate(real currentTime,
                                     real timeIncrement,
                                     size_t ltsFace,
                                     uint32_t pointIndex) {
  const real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
      currentTime - onsetTime[ltsFace][pointIndex], timeIncrement, riseTime[ltsFace][pointIndex]);
  return smoothStepIncrement / timeIncrement;
}

template <typename Cfg>
void DeltaSTF<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  onsetTime = layerData.var<LTSImposedSlipRatesDelta::OnsetTime>(Cfg());
}

template <typename Cfg>
Real<Cfg> DeltaSTF<Cfg>::evaluate(real currentTime,
                                  real timeIncrement,
                                  size_t ltsFace,
                                  uint32_t pointIndex) {
  // Currently, the delta pulse is normalized in time equivalent to FL33 and FL34
  return deltaPulse::deltaPulse(currentTime - onsetTime[ltsFace][pointIndex], timeIncrement);
}

#define SEISSOL_CONFIGITER(cfg) template class YoffeSTF<cfg>;
#include "ConfigInclude.h"
#define SEISSOL_CONFIGITER(cfg) template class GaussianSTF<cfg>;
#include "ConfigInclude.h"
#define SEISSOL_CONFIGITER(cfg) template class DeltaSTF<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::dr::friction_law::cpu
