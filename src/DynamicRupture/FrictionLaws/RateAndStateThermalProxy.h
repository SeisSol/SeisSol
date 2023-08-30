#ifndef SEISSOL_RATEANDSTATETHERMALPROXY_H
#define SEISSOL_RATEANDSTATETHERMALPROXY_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {

/*
 * Following https://eartharxiv.org/repository/view/5829/
 */
template <typename TPMethod>
class RateAndStateThermalProxyLaw
    : public RateAndStateBase<RateAndStateThermalProxyLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<RateAndStateThermalProxyLaw, TPMethod>::RateAndStateBase;

  /**
   * Copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSRateAndStateThermalPressurization const* const>(
          dynRup);

  this->srW = layerData.var(concreteLts->rsSrW);
  temperature = layerData.var(concreteLts->temperature);
  pressure = layerData.var(concreteLts->pressure);
  theta = layerData.var(concreteLts->theta);
  sigma = layerData.var(concreteLts->sigma);
  thetaTmpBuffer = layerData.var(concreteLts->thetaTmpBuffer);
  sigmaTmpBuffer = layerData.var(concreteLts->sigmaTmpBuffer);
  faultStrength = layerData.var(concreteLts->faultStrength);
  halfWidthShearZone = layerData.var(concreteLts->halfWidthShearZone);
  hydraulicDiffusivity = layerData.var(concreteLts->hydraulicDiffusivity);
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
    const double b = this->drParameters->rsB;
    const double f0 = this->drParameters->rsF0;
    const double sr0 = this->drParameters->rsSr0;
    const double localSrW = this->srW[face][pointIndex];
    const double localA = this->a[face][pointIndex];
    const double localSl0 = this->sl0[face][pointIndex];
    const double alphaTh = this->drParameters->thermalDiffusivity;
    const double alphaHy = this->hydraulicDiffusivity[face][pointIndex];
    const double specificHeat = this->drParameters->heatCapacity;
    const double lambda = this->drParameters->undrainedTPResponse;
    // TODO: add second state variable phi
    const auto phi = 0.0;
    const auto vStar = 0.0;
    const auto frictionCoefficient = 0.0;

    // steady-state state variable (eq 4)
    const auto steadyStateStateVariable = f0 + b * std::log(sr0 / localSlipRate);
    // steady-state friction coefficient (eq 6)
    // f_ss = a * arcsinh ( V / (2*V_0) * exp (psi_ss / a))
    // x in asinh(x) for mu calculation
    const double x = 0.5 / sr0 * std::exp(steadyStateStateVariable / localA) * localSlipRate;
    const auto steadyStateFrictionCoefficient = localA * std::asinh(x);
    // flash-heating steady-state friction coefficient (eq 7)
    const auto flashHeatingFrictionCoefficient = (steadyStateFrictionCoefficient - muW) / (1 + localSlipRate / localSrW) + muW;
    // flash-heating steady-state state variable (eq 8)
    const auto flashHeatingSteadyStateStateVariable = localA * std::log(2 * sr0 / localSlipRate * std::sinh( flashHeatingFrictionCoefficient / localA));
    // characteristic slip distance for TP (eq 15)
    const auto lStar = misc::power<2>(2 * specificHeat / (frictionCoefficient * lambda)) * misc::power<2>(std::sqrt(alphaHy) + std::sqrt(alphaTh)) / vStar;
    // TP proxy steady-state friction coefficient (eq 14)
    const auto tpProxyFrictionCoefficient = steadyStateFrictionCoefficient / std::cbrt( 1 + phi / lStar);
    // flash-heating steady-state state variable (eq 8)
    const auto tpProxySteadyStateStateVariable = localA * std::log(2 * sr0 / localSlipRate * std::sinh( tpProxyFrictionCoefficient / localA));

    // exact integration of dSV/dt DGL, assuming constant V over integration step
    const real exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
    const real localStateVariable =
        tpProxySteadyStateStateVariable * (1.0 - exp1) + exp1 * stateVarReference;
    assert((std::isfinite(localStateVariable) || pointIndex >= misc::numberOfBoundaryGaussPoints) &&
           "Inf/NaN detected");
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
    const double localA = this->a[ltsFace][pointIndex];
    // x in asinh(x) for mu calculation
    const double x = 0.5 / this->drParameters->rsSr0 * std::exp(localStateVariable / localA) *
                     localSlipRateMagnitude;
    const double result = localA * misc::asinh(x);
    assert((std::isfinite(result) || pointIndex >= misc::numberOfBoundaryGaussPoints) &&
           "Inf/NaN detected");
    return result;
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
    const double localA = this->a[ltsFace][pointIndex];
    const double c = 0.5 / this->drParameters->rsSr0 * std::exp(localStateVariable / localA);
    const double result =
        localA * c / std::sqrt(misc::power<2, double>(localSlipRateMagnitude * c) + 1.0);
    assert((std::isfinite(result) || pointIndex >= misc::numberOfBoundaryGaussPoints) &&
           "Inf/NaN detected");
    return result;
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
  real (*temperature)[misc::numPaddedPoints];
  real (*pressure)[misc::numPaddedPoints];
  real (*theta)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*sigma)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*thetaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*sigmaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints];
  real (*halfWidthShearZone)[misc::numPaddedPoints];
  real (*hydraulicDiffusivity)[misc::numPaddedPoints];
  real (*faultStrength)[misc::numPaddedPoints];
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_RATEANDSTATETHERMALPROXY_H
