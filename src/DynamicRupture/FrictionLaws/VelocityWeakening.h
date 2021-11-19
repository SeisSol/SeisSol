#ifndef SEISSOL_VELOCITYWEAKENING_H
#define SEISSOL_VELOCITYWEAKENING_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {
/**
 * This class was not tested and compared to the Fortran FL4. Since FL4 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
template <class Derived>
class VelocityWeakening : public RateAndStateBase<VelocityWeakening<Derived>> {
  public:
  // Attributes
  real (*stateVar)[numPaddedPoints]{};
  real (*sl0)[numPaddedPoints]{};
  real (*a)[numPaddedPoints]{};

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);
    stateVar = layerData.var(concreteLts->stateVariable);
    sl0 = layerData.var(concreteLts->rs_sl0);
    a = layerData.var(concreteLts->rs_a);
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_VELOCITYWEAKENING_H
