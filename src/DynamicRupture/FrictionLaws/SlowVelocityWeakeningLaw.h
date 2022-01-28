#ifndef SEISSOL_SLOWVELOCITYWEAKENINGLAW_H
#define SEISSOL_SLOWVELOCITYWEAKENINGLAW_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {
/**
 * This class was not tested and compared to the Fortran FL4.
 */
template <class Derived>
class SlowVelocityWeakeningLaw : public RateAndStateBase<SlowVelocityWeakeningLaw<Derived>> {
  public:
  using RateAndStateBase<SlowVelocityWeakeningLaw>::RateAndStateBase;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {}

  real updateStateVariable(int pointIndex,
                           unsigned int face,
                           real stateVarReference,
                           real timeIncrement,
                           real localSlipRate) {
    return static_cast<Derived*>(this)->updateStateVariable(
        pointIndex, face, stateVarReference, timeIncrement, localSlipRate);
  }

  /**
   * Computes the friction coefficient from the state variable and slip rate
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Theta
   * / L)}{a} \right)\right).\f]
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Theta \f$
   * @return \f$ \mu \f$
   */
  real updateMu(unsigned int ltsFace,
                unsigned int pointIndex,
                real localSlipRateMagnitude,
                real localStateVariable) {
    real localA = this->a[ltsFace][pointIndex];
    real localSl0 = this->sl0[ltsFace][pointIndex];
    real log1 = std::log(this->drParameters.rsSr0 * localStateVariable / localSl0);
    // x in asinh(x) for mu calculation
    real x = 0.5 * (localSlipRateMagnitude / this->drParameters.rsSr0) *
             std::exp((this->drParameters.rsF0 + this->drParameters.rsB * log1) / localA);
    return localA * misc::asinh(x);
  }

  /**
   * Computes the derivative of the friction coefficient with respect to the slip rate.
   * \f[\frac{\partial}{\partial V}\mu = \frac{aC}{\sqrt{(VC)^2 +1}} \text{ with } C =
   * \frac{1}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Theta / L)}{a} \right). \f]
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Theta \f$
   * @return \f$ \mu \f$
   */
  real updateMuDerivative(unsigned int ltsFace,
                          unsigned int pointIndex,
                          real localSlipRateMagnitude,
                          real localStateVariable) {
    real localA = this->a[ltsFace][pointIndex];
    real localSl0 = this->sl0[ltsFace][pointIndex];
    real log1 = std::log(this->drParameters.rsSr0 * localStateVariable / localSl0);
    real c = (0.5 / this->drParameters.rsSr0) *
             std::exp((this->drParameters.rsF0 + this->drParameters.rsB * log1) / localA);
    return localA * c / std::sqrt(misc::power<2>(localSlipRateMagnitude * c) + 1);
  }

  real calcTMP(real localStateVariable, unsigned int ltsFace, unsigned int pointIndex) {
    return 0.5 / this->drParameters.rsSr0 *
           std::exp(
               (this->drParameters.rsF0 +
                this->drParameters.rsB * std::log(this->drParameters.rsSr0 * localStateVariable /
                                                  this->sl0[ltsFace][pointIndex])) /
               this->a[ltsFace][pointIndex]);
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SLOWVELOCITYWEAKENINGLAW_H
