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
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime,
                [[maybe_unused]] real timeIncrement,
                size_t ltsFace,
                size_t pointIndex);
};

class GaussianSTF {
  private:
  real (*onsetTime)[misc::numPaddedPoints];
  real (*riseTime)[misc::numPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_SOURCETIMEFUNCTION_H
