// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/Misc.h"
#include "RateAndState.h"

#include <cmath>

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
    this->srW_ = layerData.var<LTSRateAndStateFastVelocityWeakening::RsSrW>();
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
    const double localMuW = this->muW_[faceIndex][pointIndex];
    const double localSrW = this->srW_[faceIndex][pointIndex];
    const real localA = this->a_[faceIndex][pointIndex];
    const double localSl0 = this->sl0_[faceIndex][pointIndex];

    // low-velocity steady state friction coefficient
    const real lowVelocityFriction =
        std::max(static_cast<real>(0),
                 static_cast<real>(this->f0_[faceIndex][pointIndex] -
                                   (this->b_[faceIndex][pointIndex] - localA) *
                                       log(localSlipRate / this->drParameters_->rsSr0)));
    const real steadyStateFrictionCoefficient =
        localMuW + (lowVelocityFriction - localMuW) /
                       std::pow(1.0 + misc::power<8, double>(localSlipRate / localSrW), 1.0 / 8.0);
    // TODO: check again, if double precision is necessary here (earlier, there were cancellation
    // issues)
    const real steadyStateStateVariable =
        localA * rs::logsinh(this->drParameters_->rsSr0 / localSlipRate * 2,
                             steadyStateFrictionCoefficient / localA);

    // exact integration of dSV/dt DGL, assuming constant V over integration step

    const auto preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const real exp1v = std::exp(preexp1);
    const real exp1m = -std::expm1(preexp1);
    const real localStateVariable = steadyStateStateVariable * exp1m + exp1v * stateVarReference;
    assert((std::isfinite(localStateVariable) || pointIndex >= misc::NumBoundaryGaussPoints) &&
           "Inf/NaN detected");
    return localStateVariable;
  }

  struct MuDetails {
    std::array<real, misc::NumPaddedPoints> a{};
    std::array<real, misc::NumPaddedPoints> cLin{};
    std::array<real, misc::NumPaddedPoints> cExpLog{};
    std::array<real, misc::NumPaddedPoints> cExp{};
    std::array<real, misc::NumPaddedPoints> acLin{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const real localA = this->a_[ltsFace][pointIndex];

      const real cLin = 0.5 / this->drParameters_->rsSr0;
      const real cExpLog = localStateVariable[pointIndex] / localA;
      const real cExp = rs::computeCExp(cExpLog);
      const real acLin = localA * cLin;

      details.a[pointIndex] = localA;
      details.cLin[pointIndex] = cLin;
      details.cExpLog[pointIndex] = cExpLog;
      details.cExp[pointIndex] = cExp;
      details.acLin[pointIndex] = acLin;
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
  real updateMu(std::uint32_t pointIndex, real localSlipRateMagnitude, const MuDetails& details) {
    const real lx = details.cLin[pointIndex] * localSlipRateMagnitude;
    return details.a[pointIndex] *
           rs::arsinhexp(lx, details.cExpLog[pointIndex], details.cExp[pointIndex]);
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
  real updateMuDerivative(std::uint32_t pointIndex,
                          real localSlipRateMagnitude,
                          const MuDetails& details) {
    const real lx = details.cLin[pointIndex] * localSlipRateMagnitude;
    return details.acLin[pointIndex] *
           rs::arsinhexpDerivative(lx, details.cExpLog[pointIndex], details.cExp[pointIndex]);
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
          stateVariableBuffer[pointIndex] - this->stateVariable_[ltsFace][pointIndex];
    }
    dynamicRupture::kernel::resampleParameter resampleKrnl;
    resampleKrnl.resample = init::resample::Values;
    resampleKrnl.originalQ = deltaStateVar.data();
    resampleKrnl.resampledQ = resampledDeltaStateVar.data();
    resampleKrnl.execute();

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->stateVariable_[ltsFace][pointIndex] =
          this->stateVariable_[ltsFace][pointIndex] + resampledDeltaStateVar[pointIndex];
    }
  }

  void executeIfNotConverged(const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                             std::size_t ltsFace) const {
    [[maybe_unused]] const real tmp = 0.5 / this->drParameters_->rsSr0 *
                                      std::exp(localStateVariable[0] / this->a_[ltsFace][0]) *
                                      this->slipRateMagnitude_[ltsFace][0];
    assert(!std::isnan(tmp) && "nonConvergence RS Newton");
  }

  protected:
  real (*__restrict srW_)[misc::NumPaddedPoints];
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
