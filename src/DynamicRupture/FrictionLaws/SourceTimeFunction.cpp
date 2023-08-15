#include "SourceTimeFunction.h"

namespace seissol::dr::friction_law {
template <typename Config>
void YoffeSTF<Config>::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    RealT fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesYoffe<Config> const* const>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
  tauS = layerData.var(concreteLts->tauS);
  tauR = layerData.var(concreteLts->tauR);
}

template <typename Config>
typename Config::RealT YoffeSTF<Config>::evaluate(RealT currentTime,
                                                  [[maybe_unused]] RealT timeIncrement,
                                                  size_t ltsFace,
                                                  size_t pointIndex) {
  return regularizedYoffe::regularizedYoffe(currentTime - onsetTime[ltsFace][pointIndex],
                                            tauS[ltsFace][pointIndex],
                                            tauR[ltsFace][pointIndex]);
}

template <typename Config>
void GaussianSTF<Config>::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    RealT fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesGaussian<Config> const* const>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
  riseTime = layerData.var(concreteLts->riseTime);
}

template <typename Config>
typename Config::RealT GaussianSTF<Config>::evaluate(RealT currentTime,
                                                     RealT timeIncrement,
                                                     size_t ltsFace,
                                                     size_t pointIndex) {
  const RealT smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
      currentTime - onsetTime[ltsFace][pointIndex], timeIncrement, riseTime[ltsFace][pointIndex]);
  return smoothStepIncrement / timeIncrement;
}
} // namespace seissol::dr::friction_law
