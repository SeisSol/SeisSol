#ifndef SEISSOL_STF_H
#define SEISSOL_STF_H

#include "DynamicRupture/Misc.h"
#include "Numerical_aux/RegularizedYoffe.h"
#include "Numerical_aux/GaussianNucleationFunction.h"

namespace seissol::dr::friction_law {
class YoffeSTF {
  public:
  real (*onsetTime)[misc::numPaddedPoints];
  real (*tauS)[misc::numPaddedPoints];
  real (*tauR)[misc::numPaddedPoints];

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesYoffe*>(dynRup);
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
  public:
  real (*onsetTime)[misc::numPaddedPoints];
  real (*riseTime)[misc::numPaddedPoints];

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesGaussian*>(dynRup);
    onsetTime = layerData.var(concreteLts->onsetTime);
    riseTime = layerData.var(concreteLts->riseTime);
  }

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex) {
    return gaussianNucleationFunction::smoothStepIncrement(currentTime -
                                                               onsetTime[ltsFace][pointIndex],
                                                           timeIncrement,
                                                           riseTime[ltsFace][pointIndex]) /
           timeIncrement;
  }
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_STF_H
