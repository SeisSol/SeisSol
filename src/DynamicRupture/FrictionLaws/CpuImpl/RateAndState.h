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

#include <cmath>
#include <limits>

#ifdef __INTEL_LLVM_COMPILER
#if __INTEL_LLVM_COMPILER >= 20250000
#define SEISSOL_INTEL_SIMD_EXCEPTION
#if __INTEL_LLVM_COMPILER < 20260000
#define SEISSOL_INTEL_SIMD_EXCEPTION_STRICT
#endif
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
                             const FaultStresses<Executor::Host>& initialStress,
                             TractionResults<Executor::Host>& tractionResults,
                             std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::NumPaddedPoints>& /*strengthBuffer*/,
                             std::size_t ltsFace,
                             uint32_t timeIndex) {
    bool hasConverged = false;

    // compute initial slip rate and reference values
    auto initialVariables = static_cast<Derived*>(this)->calcInitialVariables(
        faultStresses, initialStress, stateVariableBuffer, ltsFace);
    // these three are direction dependent and are swept along with the state variable below
    auto absoluteShearStress = std::move(initialVariables.absoluteShearTraction);
    auto etaInv = std::move(initialVariables.etaInv);
    auto etaNormal = std::move(initialVariables.etaNormal);
    auto slipDirection1 = std::move(initialVariables.slipDirection1);
    auto slipDirection2 = std::move(initialVariables.slipDirection2);
    auto localSlipRate = std::move(initialVariables.localSlipRate);
    auto normalStress = std::move(initialVariables.normalStress);
    auto normalStressStick = std::move(initialVariables.normalStressStick);
    const auto stateVarReference = std::move(initialVariables.stateVarReference);
    // compute slip rates by solving non-linear system of equations
    this->updateStateVariableIterative(hasConverged,
                                       stateVarReference,
                                       localSlipRate,
                                       stateVariableBuffer,
                                       normalStress,
                                       normalStressStick,
                                       absoluteShearStress,
                                       faultStresses,
                                       initialStress,
                                       etaInv,
                                       etaNormal,
                                       slipDirection1,
                                       slipDirection2,
                                       timeIndex,
                                       ltsFace);

    // compute final thermal pressure and normalStress
    tpMethod_.calcFluidPressure(
        normalStress, this->mu_, localSlipRate, this->deltaT_[timeIndex], true, ltsFace);
    updateDirectionAndProjections(slipDirection1,
                                  slipDirection2,
                                  absoluteShearStress,
                                  etaInv,
                                  etaNormal,
                                  faultStresses,
                                  initialStress,
                                  ltsFace);
    updateNormalStress(
        normalStress, normalStressStick, faultStresses, initialStress, etaNormal, ltsFace);
    // compute final slip rates and traction from average of the iterative solution and initial
    // guess
    this->calcSlipRateAndTraction(stateVarReference,
                                  localSlipRate,
                                  stateVariableBuffer,
                                  normalStress,
                                  absoluteShearStress,
                                  faultStresses,
                                  tractionResults,
                                  etaNormal,
                                  slipDirection1,
                                  slipDirection2,
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
    /// the same, before the slip rate dependent part and the clamp; the Newton solve needs it to
    /// evaluate sigma(V) itself
    std::array<real, misc::NumPaddedPoints> normalStressStick{0};
    std::array<real, misc::NumPaddedPoints> stateVarReference{0};
    std::array<real, misc::NumPaddedPoints> etaInv{0};
    std::array<real, misc::NumPaddedPoints> etaNormal{0};
    std::array<real, misc::NumPaddedPoints> slipDirection1{0};
    std::array<real, misc::NumPaddedPoints> slipDirection2{0};
  };

  /*
   * Compute shear stress magnitude, localSlipRate, effective normal stress, reference state
   * variable. Also sets slipRateMagnitude member to reference value.
   */
  InitialVariables
      calcInitialVariables(const FaultStresses<Executor::Host>& faultStresses,
                           const FaultStresses<Executor::Host>& initialStress,
                           const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                           std::size_t ltsFace) {
    // Careful, the state variable must always be corrected using stateVarZero and not
    // localStateVariable!
    std::array<real, misc::NumPaddedPoints> stateVarReference{};
    std::copy(localStateVariable.begin(), localStateVariable.end(), stateVarReference.begin());

    std::array<real, misc::NumPaddedPoints> absoluteTraction{};
    std::array<real, misc::NumPaddedPoints> normalStress{};
    std::array<real, misc::NumPaddedPoints> normalStressStick{};
    std::array<real, misc::NumPaddedPoints> temporarySlipRate{};
    std::array<real, misc::NumPaddedPoints> etaInv{};
    std::array<real, misc::NumPaddedPoints> etaNormal{};
    std::array<real, misc::NumPaddedPoints> slipDirection1{};
    std::array<real, misc::NumPaddedPoints> slipDirection2{};

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 =
          initialStress.traction1[pointIndex] + faultStresses.traction1[pointIndex];
      const real totalTraction2 =
          initialStress.traction2[pointIndex] + faultStresses.traction2[pointIndex];
      absoluteTraction[pointIndex] = misc::magnitude(totalTraction1, totalTraction2);

      const auto [_, invEta] = common::projectEta(this->impAndEta_[ltsFace],
                                                  this->impedanceMatrices_[ltsFace],
                                                  totalTraction1,
                                                  totalTraction2,
                                                  absoluteTraction[pointIndex]);

      etaInv[pointIndex] = invEta;

      etaNormal[pointIndex] = common::projectEtaNormal(this->impAndEta_[ltsFace],
                                                       this->impedanceMatrices_[ltsFace],
                                                       totalTraction1,
                                                       totalTraction2,
                                                       absoluteTraction[pointIndex]);

      // initial slip direction: the trial traction. For isotropy this stays exact.
      const real invAbsolute = (absoluteTraction[pointIndex] > 0)
                                   ? static_cast<real>(1.0) / absoluteTraction[pointIndex]
                                   : static_cast<real>(0.0);
      slipDirection1[pointIndex] = totalTraction1 * invAbsolute;
      slipDirection2[pointIndex] = totalTraction2 * invAbsolute;

      // The following process is adapted from that described by Kaneko et al. (2008)
      this->slipRateMagnitude_[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1_[ltsFace][pointIndex], this->slipRate2_[ltsFace][pointIndex]);
      this->slipRateMagnitude_[ltsFace][pointIndex] =
          std::max(rs::almostZero(), this->slipRateMagnitude_[ltsFace][pointIndex]);
      temporarySlipRate[pointIndex] = this->slipRateMagnitude_[ltsFace][pointIndex];
    } // End of pointIndex-loop

    // after the loop: updateNormalStress reads slipRateMagnitude_, which is only set above
    updateNormalStress(
        normalStress, normalStressStick, faultStresses, initialStress, etaNormal, ltsFace);

    return {absoluteTraction,
            temporarySlipRate,
            normalStress,
            normalStressStick,
            stateVarReference,
            etaInv,
            etaNormal,
            slipDirection1,
            slipDirection2};
  }

  /**
   * Anisotropy: sweeps the slip direction, and with it every direction dependent projection of the
   * impedance. A no-op for isotropy, where the slip is exactly parallel to the trial traction.
   *
   * The exact condition is tau0 = (S I + V eta_ss) n with |n| = 1, see common::updateSlipDirection.
   * The strength belonging to the current slip rate is recovered as S = n^T tau0 - V * eta_proj,
   * so no additional friction law evaluation is needed here. Riding along with the existing outer
   * fixed-point loop makes the sweep essentially free.
   */
  void updateDirectionAndProjections(
      [[maybe_unused]] std::array<real, misc::NumPaddedPoints>& slipDirection1,
      [[maybe_unused]] std::array<real, misc::NumPaddedPoints>& slipDirection2,
      [[maybe_unused]] std::array<real, misc::NumPaddedPoints>& absoluteShearTraction,
      [[maybe_unused]] std::array<real, misc::NumPaddedPoints>& etaInv,
      [[maybe_unused]] std::array<real, misc::NumPaddedPoints>& etaNormal,
      [[maybe_unused]] const FaultStresses<Executor::Host>& faultStresses,
      [[maybe_unused]] const FaultStresses<Executor::Host>& initialStress,
      [[maybe_unused]] std::size_t ltsFace) {
    if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
#pragma omp simd
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        const real totalTraction1 =
            initialStress.traction1[pointIndex] + faultStresses.traction1[pointIndex];
        const real totalTraction2 =
            initialStress.traction2[pointIndex] + faultStresses.traction2[pointIndex];
        const real trialMagnitude = misc::magnitude(totalTraction1, totalTraction2);
        const real slipRate = this->slipRateMagnitude_[ltsFace][pointIndex];

        const auto [etaProj, unusedInv] = common::projectEta(this->impAndEta_[ltsFace],
                                                             this->impedanceMatrices_[ltsFace],
                                                             slipDirection1[pointIndex],
                                                             slipDirection2[pointIndex],
                                                             static_cast<real>(1.0));

        const real strength = absoluteShearTraction[pointIndex] - slipRate * etaProj;

        const auto [n1, n2] = common::updateSlipDirection(this->impAndEta_[ltsFace],
                                                          this->impedanceMatrices_[ltsFace],
                                                          strength,
                                                          slipRate,
                                                          totalTraction1,
                                                          totalTraction2,
                                                          trialMagnitude);
        slipDirection1[pointIndex] = n1;
        slipDirection2[pointIndex] = n2;

        absoluteShearTraction[pointIndex] = n1 * totalTraction1 + n2 * totalTraction2;

        const auto [etaUnused, invEta] = common::projectEta(this->impAndEta_[ltsFace],
                                                            this->impedanceMatrices_[ltsFace],
                                                            n1,
                                                            n2,
                                                            static_cast<real>(1.0));
        etaInv[pointIndex] = invEta;
        etaNormal[pointIndex] = common::projectEtaNormal(this->impAndEta_[ltsFace],
                                                         this->impedanceMatrices_[ltsFace],
                                                         n1,
                                                         n2,
                                                         static_cast<real>(1.0));
      }
    }
  }

  void
      updateStateVariableIterative(bool& hasConverged,
                                   const std::array<real, misc::NumPaddedPoints>& stateVarReference,
                                   std::array<real, misc::NumPaddedPoints>& localSlipRate,
                                   std::array<real, misc::NumPaddedPoints>& localStateVariable,
                                   std::array<real, misc::NumPaddedPoints>& normalStress,
                                   std::array<real, misc::NumPaddedPoints>& normalStressStick,
                                   std::array<real, misc::NumPaddedPoints>& absoluteShearStress,
                                   const FaultStresses<Executor::Host>& faultStresses,
                                   const FaultStresses<Executor::Host>& initialStress,
                                   std::array<real, misc::NumPaddedPoints>& etaInv,
                                   std::array<real, misc::NumPaddedPoints>& etaNormal,
                                   std::array<real, misc::NumPaddedPoints>& slipDirection1,
                                   std::array<real, misc::NumPaddedPoints>& slipDirection2,
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

      updateDirectionAndProjections(slipDirection1,
                                    slipDirection2,
                                    absoluteShearStress,
                                    etaInv,
                                    etaNormal,
                                    faultStresses,
                                    initialStress,
                                    ltsFace);
      updateNormalStress(
          normalStress, normalStressStick, faultStresses, initialStress, etaNormal, ltsFace);

      // solve for new slip rate
      hasConverged = this->invertSlipRateIterative(ltsFace,
                                                   localStateVariable,
                                                   normalStress,
                                                   normalStressStick,
                                                   etaNormal,
                                                   absoluteShearStress,
                                                   etaInv,
                                                   testSlipRate);

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

        // Relative, like the criterion of the inner solve: the slip rates this loop walks
        // through span twenty decades, and a fixed step in m/s is a different demand at every
        // one of them -- no requirement at all where a point creeps, and more precision than
        // the inner solve resolves where it slips fast. The tolerance is raised to a few ulp,
        // since no nonzero relative step is smaller.
        const auto stateTolerance =
            std::max(static_cast<real>(this->drParameters_.rsStateTolerance),
                     static_cast<real>(4.0) * std::numeric_limits<real>::epsilon());
        const auto pointConverged =
            std::abs(testSlipRate[pointIndex] - this->slipRateMagnitude_[ltsFace][pointIndex]) <=
            stateTolerance * std::abs(testSlipRate[pointIndex]);

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
                               const std::array<real, misc::NumPaddedPoints>& absoluteTraction,
                               const FaultStresses<Executor::Host>& faultStresses,
                               TractionResults<Executor::Host>& tractionResults,
                               const std::array<real, misc::NumPaddedPoints>& etaNormal,
                               const std::array<real, misc::NumPaddedPoints>& slipDirection1,
                               const std::array<real, misc::NumPaddedPoints>& slipDirection2,
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

      // the direction along which the slip rate is decomposed; scaled such that dividing by
      // `divisor` yields the unit slip direction. For isotropy slipDirection is the (normalized)
      // trial traction and absoluteTraction is its magnitude.
      const real dirTraction1 = slipDirection1[pointIndex] * absoluteTraction[pointIndex];
      const real dirTraction2 = slipDirection2[pointIndex] * absoluteTraction[pointIndex];

      const auto [eta, _] = common::projectEta(this->impAndEta_[ltsFace],
                                               this->impedanceMatrices_[ltsFace],
                                               slipDirection1[pointIndex],
                                               slipDirection2[pointIndex],
                                               static_cast<real>(1.0));

      const auto divisor = strength + eta * this->slipRateMagnitude_[ltsFace][pointIndex];
      this->slipRate1_[ltsFace][pointIndex] =
          this->slipRateMagnitude_[ltsFace][pointIndex] * dirTraction1 / divisor;
      this->slipRate2_[ltsFace][pointIndex] =
          this->slipRateMagnitude_[ltsFace][pointIndex] * dirTraction2 / divisor;

      const auto [tU1, tU2] = common::matmulEta(this->impAndEta_[ltsFace],
                                                this->impedanceMatrices_[ltsFace],
                                                this->slipRate1_[ltsFace][pointIndex],
                                                this->slipRate2_[ltsFace][pointIndex]);

      // calculate traction
      // note that the normal stress written here is the *dynamic* normal traction, i.e. in the
      // same space as faultStresses/qInterpolated -- not the effective normal stress used for the
      // friction strength above (which additionally carries the initial stress, the initial
      // pressure, thermal pressurization and the min(., 0) clamp).
      tractionResults.normalStress[pointIndex] =
          faultStresses.normalStress[pointIndex] -
          this->slipRateMagnitude_[ltsFace][pointIndex] * etaNormal[pointIndex];
      tractionResults.traction1[pointIndex] = faultStresses.traction1[pointIndex] - tU1;
      tractionResults.traction2[pointIndex] = faultStresses.traction2[pointIndex] - tU2;
      this->traction1_[ltsFace][pointIndex] = tractionResults.traction1[pointIndex];
      this->traction2_[ltsFace][pointIndex] = tractionResults.traction2[pointIndex];

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
   * Solve for new slip rate (\f$\hat{s}\f$) with a bracketed Newton (rtsafe). \f$\hat{s}\f$ solves
   * \f[g := -\frac{1}{\eta}(|\sigma(\hat{s})|\,\mu(\hat{s}) - \Theta) - \hat{s} = 0,\f]
   * c.f. Carsten Uphoff's dissertation eq. (4.57), with \f$\eta\f$ projected onto the fault plane
   * per point and \f$\sigma\f$ following the slip rate through the anisotropic normal coupling.
   * The root is bracketed in closed form, and the bracket survives that coupling because it only
   * needs \f$\mu(0)=0\f$, \f$\mu\ge0\f$ and \f$|\sigma|\ge0\f$ (endpoints are NOT evaluated):
   *   g(0+)           = +invEta * Theta            > 0
   *   g(Theta*invEta) = -invEta * |sigma| * mu    <= 0.
   * Without the coupling \f$g' < -1\f$ and the root is unique. The coupling can weaken \f$g'\f$
   * (cf. the derivative below); bisection converges to a root inside the bracket either way, so
   * the solver does not rest on uniqueness.
   * We take Newton while it stays in the bracket and outruns bisection, else bisect. The bracket
   * is non-increasing and loses half of its decades on every fallback, so the iterate settles and
   * termination is relative in SLIP-RATE space (|dV| < xacc * V). Two floors keep that test
   * reachable in finite precision: xacc is clamped to a few ulp, and an iterate whose residual has
   * sunk into the rounding noise of its own evaluation counts as converged, because no further step
   * can be told from noise.
   */
  bool invertSlipRateIterative(std::size_t ltsFace,
                               const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                               const std::array<real, misc::NumPaddedPoints>& normalStress,
                               const std::array<real, misc::NumPaddedPoints>& normalStressStick,
                               const std::array<real, misc::NumPaddedPoints>& etaNormal,
                               const std::array<real, misc::NumPaddedPoints>& absoluteShearStress,
                               const std::array<real, misc::NumPaddedPoints>& invEta,
                               std::array<real, misc::NumPaddedPoints>& slipRateTest) {

    real muF[misc::NumPaddedPoints]{};
    real g[misc::NumPaddedPoints]{};
    real xLow[misc::NumPaddedPoints]{};
    real xHigh[misc::NumPaddedPoints]{};
    real dxOld[misc::NumPaddedPoints]{};        // previous step, for the "outrun bisection" test
    real gNoise[misc::NumPaddedPoints]{};       // rounding noise of the residual, per point
    int32_t converged[misc::NumPaddedPoints]{}; // int not bool: keeps ICX SIMD happy (cf. below)

    // Number of roundings that enter one residual evaluation; used to size both floors below.
    constexpr real NoiseFactor = 4;
    constexpr real Eps = std::numeric_limits<real>::epsilon();

    // rsSlipRateTolerance is a RELATIVE step tolerance. A nonzero step is at least one ulp of the
    // iterate, so a tolerance below Eps cannot be met at all and would leave the exact-fixed-point
    // guard as the only way out. Clamp it to a few ulp.
    const real xacc = std::max(this->drParameters_.rsSlipRateTolerance, NoiseFactor * Eps);

    const auto details = static_cast<Derived*>(this)->getMuDetails(ltsFace, localStateVariable);

    // closed-form bracket + warm start (clamped previous-step V); no endpoint evaluations
#ifndef SEISSOL_INTEL_SIMD_EXCEPTION_STRICT
#pragma omp simd
#endif
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      const real lo = rs::almostZero();
      const real hi = std::max(lo, absoluteShearStress[pointIndex] * invEta[pointIndex]);
      // A point that carries no normal stress at the free-slip limit has its root exactly there:
      // |sigma| vanishes, g is the line Theta * invEta - V, and g(hi) = 0. That is worth taking
      // directly, because it is the one root rtsafe cannot approach: a root on the bracket
      // boundary leaves the Newton step the same size as the previous one, so the guard falls
      // back to bisection on every iteration and the solve spends its whole budget halving.
      const bool openAtLimit =
          effectiveNormalStress(normalStress, normalStressStick, etaNormal, hi, pointIndex) ==
          static_cast<real>(0.0);
      xLow[pointIndex] = lo;
      xHigh[pointIndex] = hi;
      slipRateTest[pointIndex] =
          openAtLimit ? hi
                      : std::min(std::max(this->slipRateMagnitude_[ltsFace][pointIndex], lo), hi);
      dxOld[pointIndex] = hi - lo;
      converged[pointIndex] = openAtLimit ? 1 : 0;
    }

    bool allConverged = false;
    for (uint32_t i = 0; i < this->drParameters_.rsMaxNumberSlipRateUpdates; i++) {

      // residual + bracket maintenance. Transcendental mu() eval stays vectorized and unmasked;
      // converged lanes' brackets are simply never read again, so no mask is needed here.
      // >>> precision knob: evaluate muF/g (and dG below) in double -- promote |sigma| and Theta
      //     (lossless float->double) -- to drop the noise floor AND make the sign(g) test exact.
#ifndef SEISSOL_INTEL_SIMD_EXCEPTION_STRICT
#pragma omp simd
#endif
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        const real x = slipRateTest[pointIndex];
        muF[pointIndex] = static_cast<Derived*>(this)->updateMu(pointIndex, x, details);
        // sigma follows the trial slip rate, so it is evaluated at x rather than taken frozen:
        // that moves the normal coupling out of the outer fixed point and into this Newton.
        const real sigma =
            effectiveNormalStress(normalStress, normalStressStick, etaNormal, x, pointIndex);
        g[pointIndex] = -invEta[pointIndex] *
                            (std::abs(sigma) * muF[pointIndex] - absoluteShearStress[pointIndex]) -
                        x;
        // |sigma| * mu and tau cancel at the root, so the rounding error of g does not shrink
        // with the iterate: it stays at Eps times the magnitude of the two cancelling terms. Below
        // that level the sign of g -- and with it the bracket update -- carries no information.
        gNoise[pointIndex] = NoiseFactor * Eps * invEta[pointIndex] *
                             (std::abs(sigma) * muF[pointIndex] + absoluteShearStress[pointIndex]);
        const bool gPos = g[pointIndex] > static_cast<real>(0); // g decreasing: g>0 => root above x
        xLow[pointIndex] = gPos ? x : xLow[pointIndex];
        xHigh[pointIndex] = !gPos ? x : xHigh[pointIndex];
      }

      // Newton-or-bisect select (contains the division: guarded by the *non-strict* ICX macro).
#ifndef SEISSOL_INTEL_SIMD_EXCEPTION
#pragma omp simd
#endif
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        const real x = slipRateTest[pointIndex];
        const real dMuF = static_cast<Derived*>(this)->updateMuDerivative(pointIndex, x, details);
        const real sigma =
            effectiveNormalStress(normalStress, normalStressStick, etaNormal, x, pointIndex);

        // |sigma| = -sigma while the fault is closed, and sigma follows the slip rate through the
        // anisotropic normal coupling, so d|sigma|/dV = etaNormal there.
        real dAbsSigma{};
        if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
          dAbsSigma =
              (sigma < static_cast<real>(0.0)) ? etaNormal[pointIndex] : static_cast<real>(0.0);
        } else {
          dAbsSigma = static_cast<real>(0.0);
        }
        const real dGFrozen =
            -invEta[pointIndex] * (std::abs(sigma) * dMuF) - static_cast<real>(1.0);
        const real dGCoupled =
            -invEta[pointIndex] * (std::abs(sigma) * dMuF + dAbsSigma * muF[pointIndex]) -
            static_cast<real>(1.0);
        // A fault that loses normal stress as it slips (etaNormal < 0) is the only case in which
        // the coupling can weaken g. It stays strictly decreasing as long as
        //   |etaNormal| * mu < eta_proj + |sigma| * mu' ,
        // which is NOT a comfortable margin: scanning dGCoupled over the whole bracket, a locked
        // point at |etaNormal| = eta_proj already peaks at dGCoupled = -0.06, and turns positive
        // slightly above that ratio; a slipping or nearly open point holds out to roughly three
        // times eta_proj. Positive definiteness of the 3x3 impedance bounds the ratio by
        // sqrt(eta_nn / eta_proj), so a strongly anisotropic material lands close to the edge.
        // Past it g has several roots in the bracket and the solve returns one of them -- still a
        // solution of the discrete problem, but not necessarily the intended branch. This is why
        // the bracket, rather than the derivative, carries the robustness here.
        // For etaNormal > 0 the coupled derivative is the better conditioned of the two. dGFrozen
        // is negative by construction, which keeps the Newton step pointing into the bracket even
        // when the coupled one would not.
        const real dG = (dGCoupled < static_cast<real>(0.0)) ? dGCoupled : dGFrozen;

        // Which variable the residual is closer to linear in is written in dG itself: it carries
        // the strength term and the -1 of the radiation damping. Where the strength term wins,
        // mu is in its logarithmic branch and g is nearly linear in log V; where the damping
        // wins, g is nearly linear in V. The two are equal at dG = -2.
        //
        // The log step is Newton on h(u) = g(exp(u)), h'(u) = dG * V, hence multiplicative. The
        // clamp keeps exp() inside the range before the bracket test gets to reject the step, and
        // a multiplicative step cannot leave the positive axis.
        real xNewton{};
        if (dG < static_cast<real>(-2.0)) {
          const real du = std::min(std::max(-g[pointIndex] / (dG * x), static_cast<real>(-60.0)),
                                   static_cast<real>(60.0));
          xNewton = x * std::exp(du);
        } else {
          xNewton = x - g[pointIndex] / dG;
        }
        // Bisect geometrically. The bracket spans the whole admissible range of slip rates,
        // from almostZero() up to the free-slip limit tau/eta_s, so its arithmetic midpoint sits
        // many orders of magnitude above the root of a locked or creeping point, and a fallback
        // would then need one halving per factor of two to walk back down. The geometric midpoint
        // halves the number of decades instead, which is the scale the root lives on. The two
        // square roots keep the product from underflowing for the smallest brackets.
        const real xBisect = std::sqrt(xLow[pointIndex]) * std::sqrt(xHigh[pointIndex]);
        const bool useBisect =
            (xNewton <= xLow[pointIndex]) || (xNewton >= xHigh[pointIndex]) ||
            (std::abs(static_cast<real>(2.0) * g[pointIndex]) > std::abs(dxOld[pointIndex] * dG));
        const real xUpdated = useBisect ? xBisect : xNewton;
        const real step = xUpdated - x;

        const bool laneConv = (std::abs(step) < xacc * std::abs(x)) ||
                              (std::abs(g[pointIndex]) <= gNoise[pointIndex]) || (xUpdated == x);
        const int32_t nowConv = (converged[pointIndex] != 0 || laneConv) ? 1 : 0;
        converged[pointIndex] = nowConv;

        // freeze converged lanes at the accepted point (so muF stays the consistent pair);
        // advance the rest
        slipRateTest[pointIndex] = (nowConv != 0) ? x : xUpdated;
        dxOld[pointIndex] = step;
      }

      allConverged =
          std::all_of(std::begin(converged), std::end(converged), [](int32_t v) { return v != 0; });
      if (allConverged) {
        break;
      }
    }

    // publish mu consistent with the accepted slip rate
    // in case of non-convergence, flag the offending lanes -- but now
    // keyed on the x-criterion, consistent with the termination test above
#ifndef SEISSOL_INTEL_SIMD_EXCEPTION_STRICT
#pragma omp simd
#endif
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->mu_[ltsFace][pointIndex] =
          static_cast<Derived*>(this)->updateMu(pointIndex, slipRateTest[pointIndex], details);
      convergenceInner_[ltsFace][pointIndex] &= (converged[pointIndex] != 0);
    }

    return allConverged;
  }

  /**
   * Effective normal stress, including the anisotropic normal/shear coupling.
   *
   * With the 3x3 impedance eta, shear slip changes the fault-normal traction:
   *   sigma(V) = sigma_stick - V * etaNormal,   etaNormal = (eta * n)_n
   * For every isotropic material etaNormal is zero and this reduces to the previous formula.
   *
   * The slip rate is taken from slipRateMagnitude_, i.e. from the previous outer fixed-point
   * iteration (or, on entry, from the previous time step). normalStressStick keeps the part that
   * does not depend on it, so that the Newton solve can follow sigma(V) itself.
   */
  void updateNormalStress(std::array<real, misc::NumPaddedPoints>& normalStress,
                          std::array<real, misc::NumPaddedPoints>& normalStressStick,
                          const FaultStresses<Executor::Host>& faultStresses,
                          const FaultStresses<Executor::Host>& initialStress,
                          const std::array<real, misc::NumPaddedPoints>& etaNormal,
                          size_t ltsFace) {
    // Todo(SW): consider poroelastic materials together with thermal pressurization
#pragma omp simd
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      normalStressStick[pointIndex] =
          faultStresses.normalStress[pointIndex] + initialStress.normalStress[pointIndex] +
          faultStresses.fluidPressure[pointIndex] + initialStress.fluidPressure[pointIndex] -
          tpMethod_.getFluidPressure(ltsFace, pointIndex);
      normalStress[pointIndex] =
          std::min(static_cast<real>(0.0),
                   normalStressStick[pointIndex] -
                       this->slipRateMagnitude_[ltsFace][pointIndex] * etaNormal[pointIndex]);
    }
  }

  /**
   * The effective normal stress a trial slip rate belongs to.
   *
   * Without the anisotropic normal coupling this is the value updateNormalStress left behind and
   * the solve is the one every other material runs. With it, sigma follows the slip rate, and
   * evaluating it inside the Newton moves the coupling out of the outer fixed point, which
   * resolves it at a linear rate, into the quadratic one.
   */
#pragma omp declare simd
  static real
      effectiveNormalStress(const std::array<real, misc::NumPaddedPoints>& normalStress,
                            const std::array<real, misc::NumPaddedPoints>& normalStressStick,
                            const std::array<real, misc::NumPaddedPoints>& etaNormal,
                            real slipRate,
                            std::uint32_t pointIndex) {
    if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
      return std::min(static_cast<real>(0.0),
                      normalStressStick[pointIndex] - slipRate * etaNormal[pointIndex]);
    } else {
      return normalStress[pointIndex];
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
