#ifndef SEISSOL_SOURCETIMEFUNCTION_H
#define SEISSOL_SOURCETIMEFUNCTION_H

#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Numerical_aux/GaussianNucleationFunction.h"
#include "Numerical_aux/RegularizedYoffe.h"

namespace seissol::dr::friction_law {
class YoffeSTF {
  private:
  real (*onsetTime)[misc::numPaddedPoints];
  real (*tauS)[misc::numPaddedPoints];
  real (*tauR)[misc::numPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          seissol::initializer::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {
    auto* concreteLts =
        dynamic_cast<seissol::initializer::LTSImposedSlipRatesYoffe const* const>(dynRup);
    onsetTime = layerData.var(concreteLts->onsetTime);
    tauS = layerData.var(concreteLts->tauS);
    tauR = layerData.var(concreteLts->tauR);
  }
  real evaluate(real currentTime,
                [[maybe_unused]] real timeIncrement,
                size_t ltsFace,
                size_t pointIndex) {
    return regularizedYoffe::regularizedYoffe(currentTime - onsetTime[ltsFace][pointIndex],
                                              tauS[ltsFace][pointIndex],
                                              tauR[ltsFace][pointIndex]);
  }
};

class GaussianSTF {
  private:
  real (*onsetTime)[misc::numPaddedPoints];
  real (*riseTime)[misc::numPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          seissol::initializer::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {
    auto* concreteLts =
        dynamic_cast<seissol::initializer::LTSImposedSlipRatesGaussian const* const>(dynRup);
    onsetTime = layerData.var(concreteLts->onsetTime);
    riseTime = layerData.var(concreteLts->riseTime);
  }

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex) {
    const real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
        currentTime - onsetTime[ltsFace][pointIndex], timeIncrement, riseTime[ltsFace][pointIndex]);
    return smoothStepIncrement / timeIncrement;
  }
};
} // namespace seissol::dr::friction_law
#endif // SEISSOL_SOURCETIMEFUNCTION_H
