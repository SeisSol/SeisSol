#ifndef SEISSOL_SOURCETIMEFUNCTION_H
#define SEISSOL_SOURCETIMEFUNCTION_H

#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Numerical_aux/GaussianNucleationFunction.h"
#include "Numerical_aux/RegularizedYoffe.h"

namespace seissol::dr::friction_law {
template <typename Config>
class YoffeSTF {
  public:
  using RealT = typename Config::RealT;

  private:
  RealT (*onsetTime)[misc::numPaddedPoints<Config>];
  RealT (*tauS)[misc::numPaddedPoints<Config>];
  RealT (*tauR)[misc::numPaddedPoints<Config>];

  public:
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime);

  RealT evaluate(RealT currentTime,
                 [[maybe_unused]] RealT timeIncrement,
                 size_t ltsFace,
                 size_t pointIndex);
};

template <typename Config>
class GaussianSTF {
  public:
  using RealT = typename Config::RealT;

  private:
  RealT (*onsetTime)[misc::numPaddedPoints<Config>];
  RealT (*riseTime)[misc::numPaddedPoints<Config>];

  public:
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime);

  RealT evaluate(RealT currentTime, RealT timeIncrement, size_t ltsFace, size_t pointIndex);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_SOURCETIMEFUNCTION_H
