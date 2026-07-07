// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_RATEANDSTATE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_RATEANDSTATE_H_

#include "BaseFrictionLaw.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"
#include "Memory/Descriptor/DynamicRupture.h"

#ifdef __INTEL_LLVM_COMPILER
#if __INTEL_LLVM_COMPILER >= 20250000
#define SEISSOL_INTEL_SIMD_EXCEPTION
#endif
#endif // __INTEL_LLVM_COMPILER

namespace seissol::dr::friction_law::cpu {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived, class TPMethod>
class RateAndStateBase : public BaseFrictionLaw<RateAndStateBase<Derived, TPMethod>> {
  public:
  explicit RateAndStateBase(const FrictionLawParameters& drParameters)
      : BaseFrictionLaw<RateAndStateBase<Derived, TPMethod>>::BaseFrictionLaw(drParameters),
        tpMethod_(TPMethod(drParameters)) {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  void updateFrictionAndSlip(const FaultStresses<Executor::Host>& faultStresses,
                             TractionResults<Executor::Host>& tractionResults,
                             std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::NumPaddedPoints>& /*strengthBuffer*/,
                             std::size_t ltsFace,
                             uint32_t timeIndex) {
    bool hasConverged = false;

    // compute initial slip rate and reference values
    auto initialVariables = static_cast<Derived*>(this)->calcInitialVariables(
        faultStresses, stateVariableBuffer, timeIndex, ltsFace);
    const std::array<real, misc::NumPaddedPoints> absoluteShearStress =
        std::move(initialVariables.absoluteShearTraction);
    std::array<real, misc::NumPaddedPoints> localSlipRate =
        std::move(initialVariables.localSlipRate);
    std::array<real, misc::NumPaddedPoints> normalStress = std::move(initialVariables.normalStress);
    const std::array<real, misc::NumPaddedPoints> stateVarReference =
        std::move(initialVariables.stateVarReference);
    // compute slip rates by solving non-linear system of equations
    this->updateStateVariableIterative(hasConverged,
                                       stateVarReference,
                                       localSlipRate,
                                       stateVariableBuffer,
                                       normalStress,
                                       absoluteShearStress,
                                       faultStresses,
                                       timeIndex,
                                       ltsFace);

    // compute final thermal pressure and normalStress
    tpMethod_.calcFluidPressure(
        normalStress, this->mu_, localSlipRate, this->deltaT_[timeIndex], true, ltsFace);
    updateNormalStress(normalStress, faultStresses, timeIndex, ltsFace);
    // compute final slip rates and traction from average of the iterative solution and initial
    // guess
    this->calcSlipRateAndTraction(stateVarReference,
                                  localSlipRate,
                                  stateVariableBuffer,
                                  normalStress,
                                  absoluteShearStress,
                                  faultStresses,
                                  tractionResults,
                                  timeIndex,
                                  ltsFace);
  }

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {
// copy state variable from last time step
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      stateVariableBuffer[pointIndex] = this->stateVariable_[ltsFace][pointIndex];
    }
  }

  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {
    static_cast<Derived*>(this)->resampleStateVar(stateVariableBuffer, ltsFace);
  }

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {
    a_ = layerData.var<LTSRateAndState::RsA>();
    sl0_ = layerData.var<LTSRateAndState::RsSl0>();
    f0_ = layerData.var<LTSRateAndState::RsF0>();
    muW_ = layerData.var<LTSRateAndState::RsMuW>();
    b_ = layerData.var<LTSRateAndState::RsB>();
    convergenceInner_ = layerData.var<LTSRateAndState::ConvergenceInner>();
    convergenceOuter_ = layerData.var<LTSRateAndState::ConvergenceOuter>();
    stateVariable_ = layerData.var<LTSRateAndState::StateVariable>();
    static_cast<Derived*>(this)->copyStorageToLocal(layerData);
    tpMethod_.copyStorageToLocal(layerData);
  }

  /**
   * Contains all the variables, which are to be computed initially in each timestep.
   */
  struct InitialVariables {
    std::array<real, misc::NumPaddedPoints> absoluteShearTraction{0};
    std::array<real, misc::NumPaddedPoints> localSlipRate{0};
    std::array<real, misc::NumPaddedPoints> normalStress{0};
    std::array<real, misc::NumPaddedPoints> stateVarReference{0};
  };

  /*
   * Compute shear stress magnitude, localSlipRate, effective normal stress, reference state
   * variable. Also sets slipRateMagnitude member to reference value.
   */
  InitialVariables
      calcInitialVariables(const FaultStresses<Executor::Host>& faultStresses,
                           const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                           uint32_t timeIndex,
                           std::size_t ltsFace) {
    // Careful, the state variable must always be corrected using stateVarZero and not
    // localStateVariable!
    std::array<real, misc::NumPaddedPoints> stateVarReference{};
    std::copy(localStateVariable.begin(), localStateVariable.end(), stateVarReference.begin());

    std::array<real, misc::NumPaddedPoints> absoluteTraction{};
    std::array<real, misc::NumPaddedPoints> normalStress{};
    std::array<real, misc::NumPaddedPoints> temporarySlipRate{};

    updateNormalStress(normalStress, faultStresses, timeIndex, ltsFace);
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 = this->initialStressInFaultCS_[ltsFace][3][pointIndex] +
                                  faultStresses.traction1[timeIndex][pointIndex];
      const real totalTraction2 = this->initialStressInFaultCS_[ltsFace][5][pointIndex] +
                                  faultStresses.traction2[timeIndex][pointIndex];
      absoluteTraction[pointIndex] = misc::magnitude(totalTraction1, totalTraction2);

      // The following process is adapted from that described by Kaneko et al. (2008)
      this->slipRateMagnitude_[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1_[ltsFace][pointIndex], this->slipRate2_[ltsFace][pointIndex]);
      this->slipRateMagnitude_[ltsFace][pointIndex] =
          std::max(rs::almostZero(), this->slipRateMagnitude_[ltsFace][pointIndex]);
      temporarySlipRate[pointIndex] = this->slipRateMagnitude_[ltsFace][pointIndex];
    } // End of pointIndex-loop
    return {absoluteTraction, temporarySlipRate, normalStress, stateVarReference};
  }

  void updateStateVariableIterative(
      bool& hasConverged,
      const std::array<real, misc::NumPaddedPoints>& stateVarReference,
      std::array<real, misc::NumPaddedPoints>& localSlipRate,
      std::array<real, misc::NumPaddedPoints>& localStateVariable,
      std::array<real, misc::NumPaddedPoints>& normalStress,
      const std::array<real, misc::NumPaddedPoints>& absoluteShearStress,
      const FaultStresses<Executor::Host>& faultStresses,
      uint32_t timeIndex,
      std::size_t ltsFace) {
    std::array<real, misc::NumPaddedPoints> testSlipRate{};
    std::array<bool, misc::NumPaddedPoints> convergenceOuterPre{};

    // use:
    // - (inner loop) Newton-Raphson to find the fixed point slip rate with a _fixed_ state and
    // stress
    // - (outer loop) fixed-point iteration to find the fixed point slip rate with varying state and
    // stress; using the previous Newton-Raphson step (why not combine? Mainly because: thermal
    // pressurization might happen which could be a bit expensive to differentiate (though it's
    // doable; just only take all code paths that have tauV in them). Maybe the state can be
    // combined in; however that might need some more nonlinear function evaluations per step, and
    // thus could be slower). The fixed-point iteration is "regularized" by averging it with the
    // previous estimate. I.e. we compute x_(n+1) = (x_n + f(x_n)) / 2. Any fixed point we find is a
    // fixed point of f.

    // procedure source: Kaneko 2008; doi:10.1029/2007JB005154 . Section 2.3. But extended for a
    // virtually unlimited number of outer fixed point iterations.

    for (uint32_t j = 0; j < this->drParameters_.rsNumberStateVariableUpdates; j++) {
#pragma omp simd
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // fault strength using friction coefficient and fluid pressure from previous
        // timestep/iteration update state variable using sliprate from the previous time step
        localStateVariable[pointIndex] =
            static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                             ltsFace,
                                                             stateVarReference[pointIndex],
                                                             this->deltaT_[timeIndex],
                                                             localSlipRate[pointIndex]);
      }
      tpMethod_.calcFluidPressure(
          normalStress, this->mu_, localSlipRate, this->deltaT_[timeIndex], false, ltsFace);

      updateNormalStress(normalStress, faultStresses, timeIndex, ltsFace);

      // solve for new slip rate
      hasConverged = this->invertSlipRateIterative(
          ltsFace, localStateVariable, normalStress, absoluteShearStress, testSlipRate);

      // int for ICX not to fail
      int32_t converged = 1;

#pragma omp simd reduction(min : converged)
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // update local slip rate, now using V=(Vnew+Vold)/2
        // For the next SV update, use the mean slip rate between the initial guess and the one
        // found (Kaneko 2008, step 6)
        localSlipRate[pointIndex] =
            static_cast<real>(0.5) *
            (this->slipRateMagnitude_[ltsFace][pointIndex] + testSlipRate[pointIndex]);

        const auto pointConverged =
            std::abs(testSlipRate[pointIndex] - this->slipRateMagnitude_[ltsFace][pointIndex]) <
            this->drParameters_.rsStateTolerance;

        converged = std::min(pointConverged ? 1 : 0, converged);

        convergenceOuterPre[pointIndex] = pointConverged;

        // solve again for Vnew
        this->slipRateMagnitude_[ltsFace][pointIndex] = testSlipRate[pointIndex];
      } // End of pointIndex-loop

      if (converged != 0) {
        break;
      }
    }

    // update (outer) (non-)convergence
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      convergenceOuter_[ltsFace][pointIndex] &= convergenceOuterPre[pointIndex];
    }
  }

  void calcSlipRateAndTraction(const std::array<real, misc::NumPaddedPoints>& stateVarReference,
                               const std::array<real, misc::NumPaddedPoints>& localSlipRate,
                               std::array<real, misc::NumPaddedPoints>& localStateVariable,
                               const std::array<real, misc::NumPaddedPoints>& normalStress,
                               const std::array<real, misc::NumPaddedPoints>& /*absoluteTraction*/,
                               const FaultStresses<Executor::Host>& faultStresses,
                               TractionResults<Executor::Host>& tractionResults,
                               uint32_t timeIndex,
                               std::size_t ltsFace) {

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // SV from mean slip rate in tmp
      localStateVariable[pointIndex] =
          static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                           ltsFace,
                                                           stateVarReference[pointIndex],
                                                           this->deltaT_[timeIndex],
                                                           localSlipRate[pointIndex]);
    }

    const auto details = static_cast<Derived*>(this)->getMuDetails(ltsFace, localStateVariable);

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // update LocMu for next strength determination, only needed for last update
      this->mu_[ltsFace][pointIndex] = static_cast<Derived*>(this)->updateMu(
          pointIndex, this->slipRateMagnitude_[ltsFace][pointIndex], details);
      const real strength = -this->mu_[ltsFace][pointIndex] * normalStress[pointIndex];
      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 = this->initialStressInFaultCS_[ltsFace][3][pointIndex] +
                                  faultStresses.traction1[timeIndex][pointIndex];
      const real totalTraction2 = this->initialStressInFaultCS_[ltsFace][5][pointIndex] +
                                  faultStresses.traction2[timeIndex][pointIndex];

      const auto divisor =
          strength + this->impAndEta_[ltsFace].etaS * this->slipRateMagnitude_[ltsFace][pointIndex];
      this->slipRate1_[ltsFace][pointIndex] =
          this->slipRateMagnitude_[ltsFace][pointIndex] * totalTraction1 / divisor;
      this->slipRate2_[ltsFace][pointIndex] =
          this->slipRateMagnitude_[ltsFace][pointIndex] * totalTraction2 / divisor;

      // calculate traction
      tractionResults.traction1[timeIndex][pointIndex] =
          faultStresses.traction1[timeIndex][pointIndex] -
          this->impAndEta_[ltsFace].etaS * this->slipRate1_[ltsFace][pointIndex];
      tractionResults.traction2[timeIndex][pointIndex] =
          faultStresses.traction2[timeIndex][pointIndex] -
          this->impAndEta_[ltsFace].etaS * this->slipRate2_[ltsFace][pointIndex];
      this->traction1_[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
      this->traction2_[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];

      // Compute slip
      // ABS of locSlipRate removed as it would be the accumulated slip that is usually not needed
      // in the solver, see linear slip weakening
      this->accumulatedSlipMagnitude_[ltsFace][pointIndex] +=
          this->slipRateMagnitude_[ltsFace][pointIndex] * this->deltaT_[timeIndex];

      // update directional slip
      this->slip1_[ltsFace][pointIndex] +=
          this->slipRate1_[ltsFace][pointIndex] * this->deltaT_[timeIndex];
      this->slip2_[ltsFace][pointIndex] +=
          this->slipRate2_[ltsFace][pointIndex] * this->deltaT_[timeIndex];
    }
  }

  void saveDynamicStressOutput(std::size_t faceIndex, real time) {
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {

      if (this->ruptureTime_[faceIndex][pointIndex] > static_cast<real>(0.0) &&
          this->ruptureTime_[faceIndex][pointIndex] <= time &&
          this->dynStressTimePending_[faceIndex][pointIndex] &&
          this->mu_[faceIndex][pointIndex] <=
              (this->muW_[faceIndex][pointIndex] +
               static_cast<real>(0.05) *
                   (this->f0_[faceIndex][pointIndex] - this->muW_[faceIndex][pointIndex]))) {
        this->dynStressTime_[faceIndex][pointIndex] = time;
        this->dynStressTimePending_[faceIndex][pointIndex] = false;
      }
    }
  }

  /**
   * Solve for new slip rate (\f$\hat{s}\f$), applying the Newton-Raphson algorithm.
   * \f$\hat{s}\f$ has to fulfill
   * \f[g := \frac{1}{\eta_s} \cdot (\sigma_n \cdot \mu - \Theta) - \hat{s} = 0.\f] c.f. Carsten
   * Uphoff's dissertation eq. (4.57). Find root of \f$g\f$ with \f$g^\prime = \partial g / \partial
   * \hat{s}\f$: \f$\hat{s}_{i+1} = \hat{s}_i - ( g_i / g^\prime_i )\f$
   * @param ltsFace index of the faceIndex for which we invert the sliprate
   * @param localStateVariable \f$\psi\f$, needed to compute \f$\mu = f(\hat{s}, \psi)\f$
   * @param normalStress \f$\sigma_n\f$
   * @param absoluteShearStress \f$\Theta\f$
   * @param slipRateTest \f$\hat{s}\f$
   */
  bool invertSlipRateIterative(std::size_t ltsFace,
                               const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                               const std::array<real, misc::NumPaddedPoints>& normalStress,
                               const std::array<real, misc::NumPaddedPoints>& absoluteShearStress,
                               std::array<real, misc::NumPaddedPoints>& slipRateTest) {

    real muF[misc::NumPaddedPoints]{};
    real g[misc::NumPaddedPoints]{};

    const auto details = static_cast<Derived*>(this)->getMuDetails(ltsFace, localStateVariable);

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // first guess = sliprate value of the previous step
      slipRateTest[pointIndex] = this->slipRateMagnitude_[ltsFace][pointIndex];
    }

    for (uint32_t i = 0; i < this->drParameters_.rsMaxNumberSlipRateUpdates; i++) {

#ifdef SEISSOL_INTEL_SIMD_EXCEPTION
#pragma clang loop vectorize(enable)
#else
#pragma omp simd
#endif
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // calculate friction coefficient and objective function
        muF[pointIndex] =
            static_cast<Derived*>(this)->updateMu(pointIndex, slipRateTest[pointIndex], details);
        g[pointIndex] = -this->impAndEta_[ltsFace].invEtaS *
                            (std::abs(normalStress[pointIndex]) * muF[pointIndex] -
                             absoluteShearStress[pointIndex]) -
                        slipRateTest[pointIndex];
      }

      // max element of g must be smaller than newtonTolerance
      const bool hasConverged = std::all_of(std::begin(g), std::end(g), [&](auto val) {
        return std::abs(val) < this->drParameters_.rsSlipRateTolerance;
      });
      if (hasConverged) {
#pragma omp simd
        for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
          this->mu_[ltsFace][pointIndex] = muF[pointIndex];
        }
        return hasConverged;
      }

#ifdef SEISSOL_INTEL_SIMD_EXCEPTION
#pragma clang loop vectorize(enable)
#else
#pragma omp simd
#endif
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        const auto localSlipRateTest = slipRateTest[pointIndex];

        const auto dMuF =
            static_cast<Derived*>(this)->updateMuDerivative(pointIndex, localSlipRateTest, details);

        // derivative of g
        const auto dG =
            -this->impAndEta_[ltsFace].invEtaS * (std::abs(normalStress[pointIndex]) * dMuF) -
            static_cast<real>(1.0);

        // newton update
        const real tmp3 = g[pointIndex] / dG;
        slipRateTest[pointIndex] = std::max(rs::almostZero(), localSlipRateTest - tmp3);
      }
    }

    // update (non-)convergence
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      convergenceInner_[ltsFace][pointIndex] &=
          std::abs(g[pointIndex]) < this->drParameters_.rsSlipRateTolerance;
    }
    return false;
  }

  void updateNormalStress(std::array<real, misc::NumPaddedPoints>& normalStress,
                          const FaultStresses<Executor::Host>& faultStresses,
                          size_t timeIndex,
                          size_t ltsFace) {
    // Todo(SW): consider poroelastic materials together with thermal pressurisation
#pragma omp simd
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      normalStress[pointIndex] =
          std::min(static_cast<real>(0.0),
                   faultStresses.normalStress[timeIndex][pointIndex] +
                       this->initialStressInFaultCS_[ltsFace][0][pointIndex] +
                       faultStresses.fluidPressure[timeIndex][pointIndex] +
                       this->initialPressure_[ltsFace][pointIndex] -
                       tpMethod_.getFluidPressure(ltsFace, pointIndex));
    }
  }

  protected:
  // Attributes
  real (*__restrict a_)[misc::NumPaddedPoints]{};
  real (*__restrict sl0_)[misc::NumPaddedPoints]{};
  real (*__restrict stateVariable_)[misc::NumPaddedPoints]{};

  real (*__restrict f0_)[misc::NumPaddedPoints]{};
  real (*__restrict muW_)[misc::NumPaddedPoints]{};
  real (*__restrict b_)[misc::NumPaddedPoints]{};

  bool (*__restrict convergenceInner_)[misc::NumPaddedPoints]{};
  bool (*__restrict convergenceOuter_)[misc::NumPaddedPoints]{};

  TPMethod tpMethod_;
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_RATEANDSTATE_H_
