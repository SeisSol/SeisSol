#include "SourceTimeFunction.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/Tree/Layer.h"
#include "Kernels/Precision.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"
#include <cstddef>

namespace seissol::dr::friction_law {
void YoffeSTF::copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                                  const seissol::initializer::DynamicRupture* const dynRup,
                                  real fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesYoffe* const>(dynRup);
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
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesGaussian* const>(dynRup);
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
  auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSImposedSlipRatesDelta* const>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
  godunovData = layerData.var(concreteLts->godunovData);
}

real DeltaSTF::evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex) {
  // Currently, the delta pulse is normalized in time equivalent to FL33 and FL34
  const double surfaceArea = godunovData[ltsFace].doubledSurfaceArea * 0.5;
  return deltaPulse::deltaPulse(currentTime - onsetTime[ltsFace][pointIndex], timeIncrement, surfaceArea);
}

} // namespace seissol::dr::friction_law
