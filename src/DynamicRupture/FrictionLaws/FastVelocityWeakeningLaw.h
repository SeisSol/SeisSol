#ifndef SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H
#define SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {

template <typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /**
   * Copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {
    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSRateAndStateFastVelocityWeakening const* const>(
            dynRup);

    this->srW = layerData.var(concreteLts->rsSrW);
  }

  /**
   * Integrates the state variable ODE in time
   * \f[\frac{\partial \Psi}{\partial t} = - \frac{V}{L}\left(\Psi - \Psi_{ss}(V) \right)\f]
   * with steady state variable \f$\Psi_{ss}\f$.
   * Assume \f$V\f$ is constant through the time interval, then the analytic solution is:
   * \f[ \Psi(t) = \Psi_0 \exp\left( -\frac{V}{L} t \right) + \Psi_{ss} \left( 1 - \exp\left(
   * - \frac{V}{L} t\right) \right).\f]
   * @param stateVarReference \f$ \Psi_0 \f$
   * @param timeIncrement \f$ t \f$
   * @param localSlipRate \f$ V \f$
   * @return \f$ \Psi(t) \f$
   */
  #pragma omp declare simd
  real updateStateVariable(unsigned int pointIndex,
                           unsigned int face,
                           real stateVarReference,
                           real timeIncrement,
                           real localSlipRate) const {
    const double muW = this->drParameters->muW;
    const double localSrW = this->srW[face][pointIndex];
    const double localA = this->a[face][pointIndex];
    const double localSl0 = this->sl0[face][pointIndex];

    // low-velocity steady state friction coefficient
    const real lowVelocityFriction =
        this->drParameters->rsF0 -
        (this->drParameters->rsB - localA) * log(localSlipRate / this->drParameters->rsSr0);
    const real steadyStateFrictionCoefficient =
        muW + (lowVelocityFriction - muW) /
                  std::pow(1.0 + misc::power<8, double>(localSlipRate / localSrW), 1.0 / 8.0);
    // For compiling reasons we write SINH(X)=(EXP(X)-EXP(-X))/2
    // Need double precision here, otherwise this will evaluate to infinity
    const real steadyStateStateVariable =
        localA * std::log(this->drParameters->rsSr0 / localSlipRate *
                          (std::exp(steadyStateFrictionCoefficient / localA) -
                           std::exp(-steadyStateFrictionCoefficient / localA)));

    // exact integration of dSV/dt DGL, assuming constant V over integration step

    const real exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
    const real localStateVariable =
        steadyStateStateVariable * (1.0 - exp1) + exp1 * stateVarReference;
    assert(!(std::isnan(localStateVariable) && pointIndex < misc::numberOfBoundaryGaussPoints) &&
           "NaN detected");
    return localStateVariable;
  }

  /**
   * Computes the friction coefficient from the state variable and slip rate
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp
   * \left(\frac{\Psi}{a}\right)\right). \f]
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Psi \f$
   * @return \f$ \mu \f$
   */
  #pragma omp declare simd
  real updateMu(unsigned int ltsFace,
                unsigned int pointIndex,
                real localSlipRateMagnitude,
                real localStateVariable) const {
    // mu = a * arcsinh ( V / (2*V_0) * exp (psi / a))
    const real localA = this->a[ltsFace][pointIndex];
    // x in asinh(x) for mu calculation
    const real x = 0.5 / this->drParameters->rsSr0 * std::exp(localStateVariable / localA) *
                   localSlipRateMagnitude;
    return localA * misc::asinh(x);
  }

  /**
   * Computes the derivative of the friction coefficient with respect to the slip rate.
   * \f[\frac{\partial}{\partial V}\mu = \frac{aC}{\sqrt{ (VC)^2 + 1} \text{ with } C =
   * \frac{1}{2V_0} \cdot \exp \left(\frac{\Psi}{a}\right)\right).\f]
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Psi \f$
   * @return \f$ \mu \f$
   */
  #pragma omp declare simd
  real updateMuDerivative(unsigned int ltsFace,
                          unsigned int pointIndex,
                          real localSlipRateMagnitude,
                          real localStateVariable) const {
    const real localA = this->a[ltsFace][pointIndex];
    const real c = 0.5 / this->drParameters->rsSr0 * std::exp(localStateVariable / localA);
    return localA * c / std::sqrt(misc::power<2, double>(localSlipRateMagnitude * c) + 1.0);
  }

  /**
   * Resample the state variable.
   */
  void resampleStateVar(std::array<real, misc::numPaddedPoints> const& stateVariableBuffer,
                        unsigned int ltsFace) const {
    std::array<real, misc::numPaddedPoints> deltaStateVar = {0};
    std::array<real, misc::numPaddedPoints> resampledDeltaStateVar = {0};
    #pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      deltaStateVar[pointIndex] =
          stateVariableBuffer[pointIndex] - this->stateVariable[ltsFace][pointIndex];
    }
    dynamicRupture::kernel::resampleParameter resampleKrnl;
    resampleKrnl.resample = init::resample::Values;
    resampleKrnl.originalQ = deltaStateVar.data();
    resampleKrnl.resampledQ = resampledDeltaStateVar.data();
    resampleKrnl.execute();

    #pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      this->stateVariable[ltsFace][pointIndex] =
          this->stateVariable[ltsFace][pointIndex] + resampledDeltaStateVar[pointIndex];
    }
  }

  void executeIfNotConverged(std::array<real, misc::numPaddedPoints> const& localStateVariable,
                             unsigned ltsFace) const {
    [[maybe_unused]] const real tmp = 0.5 / this->drParameters->rsSr0 *
                                      exp(localStateVariable[0] / this->a[ltsFace][0]) *
                                      this->slipRateMagnitude[ltsFace][0];
    assert(!std::isnan(tmp) && "nonConvergence RS Newton");
  }

  protected:
  real (*srW)[misc::numPaddedPoints];
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_RATEANDSTATEFASTVELOCITYWEAKENING_H
