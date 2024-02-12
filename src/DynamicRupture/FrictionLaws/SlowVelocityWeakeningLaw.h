#ifndef SEISSOL_SLOWVELOCITYWEAKENINGLAW_H
#define SEISSOL_SLOWVELOCITYWEAKENINGLAW_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {
template <class Derived, class TPMethod>
class SlowVelocityWeakeningLaw
    : public RateAndStateBase<SlowVelocityWeakeningLaw<Derived, TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<SlowVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {}

// Note that we need double precision here, since single precision led to NaNs.
#pragma omp declare simd
  double updateStateVariable(int pointIndex,
                             unsigned int face,
                             double stateVarReference,
                             double timeIncrement,
                             double localSlipRate) {
    return static_cast<Derived*>(this)->updateStateVariable(
        pointIndex, face, stateVarReference, timeIncrement, localSlipRate);
  }

/**
 * Computes the friction coefficient from the state variable and slip rate
 * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Psi
 * / L)}{a} \right)\right).\f]
 * Note that we need double precision here, since single precision led to NaNs.
 * @param localSlipRateMagnitude \f$ V \f$
 * @param localStateVariable \f$ \Psi \f$
 * @return \f$ \mu \f$
 */
#pragma omp declare simd
  double updateMu(unsigned int ltsFace,
                  unsigned int pointIndex,
                  double localSlipRateMagnitude,
                  double localStateVariable) {
    const double localA = this->a[ltsFace][pointIndex];
    const double localSl0 = this->sl0[ltsFace][pointIndex];
    const double log1 = std::log(this->drParameters->rsSr0 * localStateVariable / localSl0);
    // x in asinh(x) for mu calculation
    const double x = 0.5 * (localSlipRateMagnitude / this->drParameters->rsSr0) *
                     std::exp((this->drParameters->rsF0 + this->drParameters->rsB * log1) / localA);
    return localA * misc::asinh(x);
  }

/**
 * Computes the derivative of the friction coefficient with respect to the slip rate.
 * \f[\frac{\partial}{\partial V}\mu = \frac{aC}{\sqrt{(VC)^2 +1}} \text{ with } C =
 * \frac{1}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Psi / L)}{a} \right). \f]
 * Note that we need double precision here, since single precision led to NaNs.
 * @param localSlipRateMagnitude \f$ V \f$
 * @param localStateVariable \f$ \Psi \f$
 * @return \f$ \mu \f$
 */
#pragma omp declare simd
  double updateMuDerivative(unsigned int ltsFace,
                            unsigned int pointIndex,
                            double localSlipRateMagnitude,
                            double localStateVariable) {
    const double localA = this->a[ltsFace][pointIndex];
    const double localSl0 = this->sl0[ltsFace][pointIndex];
    const double log1 = std::log(this->drParameters->rsSr0 * localStateVariable / localSl0);
    const double c = (0.5 / this->drParameters->rsSr0) *
                     std::exp((this->drParameters->rsF0 + this->drParameters->rsB * log1) / localA);
    return localA * c / std::sqrt(misc::power<2>(localSlipRateMagnitude * c) + 1);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws, we just copy the buffer into the
   * member variable.
   */
  void resampleStateVar(std::array<real, misc::numPaddedPoints> const& stateVariableBuffer,
                        unsigned int ltsFace) const {
#pragma omp simd
    for (size_t pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      this->stateVariable[ltsFace][pointIndex] = stateVariableBuffer[pointIndex];
    }
  }

  void executeIfNotConverged(std::array<real, misc::numPaddedPoints> const& localStateVariable,
                             unsigned ltsFace) {
    [[maybe_unused]] const real tmp =
        0.5 / this->drParameters->rsSr0 *
        std::exp(
            (this->drParameters->rsF0 +
             this->drParameters->rsB * std::log(this->drParameters->rsSr0 * localStateVariable[0] /
                                                this->drParameters->rsSr0)) /
            this->a[ltsFace][0]);
    assert(!std::isnan(tmp) && "nonConvergence RS Newton");
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SLOWVELOCITYWEAKENINGLAW_H
