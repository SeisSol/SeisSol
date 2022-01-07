#ifndef SEISSOL_SLOWVELOCITYWEAKENINGLAW_H
#define SEISSOL_SLOWVELOCITYWEAKENINGLAW_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {
/**
 * This class was not tested and compared to the Fortran FL4. Since FL4 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
template <class Derived>
class SlowVelocityWeakeningLaw : public RateAndStateBase<SlowVelocityWeakeningLaw<Derived>> {
  public:
  using RateAndStateBase<SlowVelocityWeakeningLaw>::RateAndStateBase;
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

  real updateStateVariable(int pointIndex,
                           unsigned int face,
                           real stateVarReference,
                           real timeIncrement,
                           real& localSlipRate) {
    return static_cast<Derived*>(this)->updateStateVariable(
        pointIndex, face, stateVarReference, timeIncrement, localSlipRate);
  }

  /**
   * Computes the friction coefficient from the state variable and slip rate
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Theta
   * / L)}{a} \right)\right).\f] \f$V\f$ is taken from the stored values.
   * @param localStateVariable \f$ \Theta \f$
   * @return \f$ \mu \f$
   */
  real updateMu(unsigned int ltsFace, unsigned int pointIndex, real localStateVariable) {
    real localA = a[ltsFace][pointIndex];
    real localSL0 = sl0[ltsFace][pointIndex];
    real log1 = std::log(this->drParameters.rs_sr0 * localStateVariable / localSL0);
    // x in asinh(x) for mu calculation
    real x = 0.5 * (this->slipRateMagnitude[ltsFace][pointIndex] / this->drParameters.rs_sr0) *
             exp((this->drParameters.rs_f0 + this->drParameters.rs_b * log1) / localA);
    return localA * std::log(x + std::sqrt(std::pow(x, 2) + 1.0));
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SLOWVELOCITYWEAKENINGLAW_H
