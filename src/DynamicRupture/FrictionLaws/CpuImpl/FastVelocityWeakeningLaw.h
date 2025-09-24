// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_

#include <DynamicRupture/Misc.h>
#include <cmath>

#include "RateAndState.h"

namespace seissol::dr::friction_law::cpu {

template <typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /**
   * Copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyStorageToLocal(DynamicRupture::Layer& layerData) {
    this->srW = layerData.var<LTSRateAndStateFastVelocityWeakening::RsSrW>();
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
  [[nodiscard]] real updateStateVariable(std::uint32_t pointIndex,
                                         std::size_t faceIndex,
                                         real stateVarReference,
                                         real timeIncrement,
                                         real localSlipRate) const {
    const double localMuW = this->muW[faceIndex][pointIndex];
    const double localSrW = this->srW[faceIndex][pointIndex];
    const double localA = this->a[faceIndex][pointIndex];
    const double localSl0 = this->sl0[faceIndex][pointIndex];

    // low-velocity steady state friction coefficient
    const real lowVelocityFriction =
        this->f0[faceIndex][pointIndex] -
        (this->b[faceIndex][pointIndex] - localA) * log(localSlipRate / this->drParameters->rsSr0);
    const real steadyStateFrictionCoefficient =
        localMuW + (lowVelocityFriction - localMuW) /
                       std::pow(1.0 + misc::power<8, double>(localSlipRate / localSrW), 1.0 / 8.0);
    // TODO: check again, if double precision is necessary here (earlier, there were cancellation
    // issues)
    const real steadyStateStateVariable =
        localA * std::log(this->drParameters->rsSr0 / localSlipRate * 2 *
                          std::sinh(steadyStateFrictionCoefficient / localA));

    // exact integration of dSV/dt DGL, assuming constant V over integration step

    const auto preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const real exp1 = std::exp(preexp1);
    const real exp1m = -std::expm1(preexp1);
    const real localStateVariable = steadyStateStateVariable * exp1m + exp1 * stateVarReference;
    assert((std::isfinite(localStateVariable) || pointIndex >= misc::NumBoundaryGaussPoints) &&
           "Inf/NaN detected");
    return localStateVariable;
  }

  struct MuDetails {
    std::array<double, misc::NumPaddedPoints> a{};
    std::array<double, misc::NumPaddedPoints> c{};
    std::array<double, misc::NumPaddedPoints> ac{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const double localA = this->a[ltsFace][pointIndex];
      const double c = 0.5 / this->drParameters->rsSr0 *
                       std::exp(static_cast<double>(localStateVariable[pointIndex]) / localA);
      assert((std::isfinite(c) || pointIndex >= misc::NumBoundaryGaussPoints) &&
             "Inf/NaN detected");
      details.a[pointIndex] = localA;
      details.c[pointIndex] = c;
      details.ac[pointIndex] = localA * c;
    }
    return details;
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
  double
      updateMu(std::uint32_t pointIndex, double localSlipRateMagnitude, const MuDetails& details) {
    const double x = details.c[pointIndex] * localSlipRateMagnitude;
    return details.a[pointIndex] * std::asinh(x);
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
  double updateMuDerivative(std::uint32_t pointIndex,
                            double localSlipRateMagnitude,
                            const MuDetails& details) {
    const double x = details.c[pointIndex] * localSlipRateMagnitude;
    return details.ac[pointIndex] / std::sqrt(x * x + 1.0);
  }

  /**
   * Resample the state variable.
   */
  void resampleStateVar(const std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                        std::size_t ltsFace) const {
    std::array<real, misc::NumPaddedPoints> deltaStateVar = {0};
    std::array<real, misc::NumPaddedPoints> resampledDeltaStateVar = {0};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      deltaStateVar[pointIndex] =
          stateVariableBuffer[pointIndex] - this->stateVariable[ltsFace][pointIndex];
    }
    dynamicRupture::kernel::resampleParameter resampleKrnl;
    resampleKrnl.resample = init::resample::Values;
    resampleKrnl.originalQ = deltaStateVar.data();
    resampleKrnl.resampledQ = resampledDeltaStateVar.data();
    resampleKrnl.execute();

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->stateVariable[ltsFace][pointIndex] =
          this->stateVariable[ltsFace][pointIndex] + resampledDeltaStateVar[pointIndex];
    }
  }

  void executeIfNotConverged(const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                             std::size_t ltsFace) const {
    [[maybe_unused]] const real tmp = 0.5 / this->drParameters->rsSr0 *
                                      std::exp(localStateVariable[0] / this->a[ltsFace][0]) *
                                      this->slipRateMagnitude[ltsFace][0];
    assert(!std::isnan(tmp) && "nonConvergence RS Newton");
  }

  protected:
  real (*__restrict srW)[misc::NumPaddedPoints];
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
