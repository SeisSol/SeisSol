#include "SourceTimeFunction.h"
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/Layer.hpp>
#include <Kernels/precision.hpp>
#include <Numerical_aux/GaussianNucleationFunction.h>
#include <Numerical_aux/RegularizedYoffe.h>
#include <cstddef>

namespace seissol::dr::friction_law {
void YoffeSTF::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                  seissol::initializer::DynamicRupture const* const dynRup,
                                  real fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<seissol::initializer::LTSImposedSlipRatesYoffe const* const>(dynRup);
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
                                     seissol::initializer::DynamicRupture const* const dynRup,
                                     real fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<seissol::initializer::LTSImposedSlipRatesGaussian const* const>(dynRup);
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
} // namespace seissol::dr::friction_law
