#ifndef SEISSOL_SOURCETIMEFUNCTION_H
#define SEISSOL_SOURCETIMEFUNCTION_H

#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"

namespace seissol::dr::friction_law {
class YoffeSTF {
  private:
  real (*onsetTime)[misc::NumPaddedPoints];
  real (*tauS)[misc::NumPaddedPoints];
  real (*tauR)[misc::NumPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime,
                [[maybe_unused]] real timeIncrement,
                size_t ltsFace,
                size_t pointIndex);
};

class GaussianSTF {
  private:
  real (*onsetTime)[misc::NumPaddedPoints];
  real (*riseTime)[misc::NumPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_SOURCETIMEFUNCTION_H
