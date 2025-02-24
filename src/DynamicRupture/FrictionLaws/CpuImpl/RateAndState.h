// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATE_H_

#include "BaseFrictionLaw.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"

namespace seissol::dr::friction_law::cpu {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived, class TPMethod>
class RateAndStateBase : public BaseFrictionLaw<RateAndStateBase<Derived, TPMethod>> {
  public:
  explicit RateAndStateBase(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionLaw<RateAndStateBase<Derived, TPMethod>>::BaseFrictionLaw(drParameters),
        tpMethod(TPMethod(drParameters)) {}

  void updateFrictionAndSlip(const FaultStresses<Executor::Host>& faultStresses,
                             TractionResults<Executor::Host>& tractionResults,
                             std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::NumPaddedPoints>& strengthBuffer,
                             unsigned ltsFace,
                             unsigned timeIndex) {
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

    // check for convergence
    if (!hasConverged) {
      static_cast<Derived*>(this)->executeIfNotConverged(stateVariableBuffer, ltsFace);
    }
    // compute final thermal pressure and normalStress
    tpMethod.calcFluidPressure(
        normalStress, this->mu, localSlipRate, this->deltaT[timeIndex], true, timeIndex, ltsFace);
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

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {
// copy state variable from last time step
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      stateVariableBuffer[pointIndex] = this->stateVariable[ltsFace][pointIndex];
    }
  }

  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {
    static_cast<Derived*>(this)->resampleStateVar(stateVariableBuffer, ltsFace);
  }

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {
    const auto* concreteLts = dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup);
    a = layerData.var(concreteLts->rsA);
    sl0 = layerData.var(concreteLts->rsSl0);
    stateVariable = layerData.var(concreteLts->stateVariable);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    tpMethod.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
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
                           unsigned int timeIndex,
                           unsigned int ltsFace) {
    // Careful, the state variable must always be corrected using stateVarZero and not
    // localStateVariable!
    std::array<real, misc::NumPaddedPoints> stateVarReference{};
    std::copy(localStateVariable.begin(), localStateVariable.end(), stateVarReference.begin());

    std::array<real, misc::NumPaddedPoints> absoluteTraction{};
    std::array<real, misc::NumPaddedPoints> normalStress{};
    std::array<real, misc::NumPaddedPoints> temporarySlipRate{};

    updateNormalStress(normalStress, faultStresses, timeIndex, ltsFace);
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                                  faultStresses.traction1[timeIndex][pointIndex];
      const real totalTraction2 = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                                  faultStresses.traction2[timeIndex][pointIndex];
      absoluteTraction[pointIndex] = misc::magnitude(totalTraction1, totalTraction2);

      // The following process is adapted from that described by Kaneko et al. (2008)
      this->slipRateMagnitude[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1[ltsFace][pointIndex], this->slipRate2[ltsFace][pointIndex]);
      this->slipRateMagnitude[ltsFace][pointIndex] =
          std::max(rs::almostZero(), this->slipRateMagnitude[ltsFace][pointIndex]);
      temporarySlipRate[pointIndex] = this->slipRateMagnitude[ltsFace][pointIndex];
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
      unsigned int timeIndex,
      unsigned int ltsFace) {
    std::array<real, misc::NumPaddedPoints> testSlipRate{0};
    for (unsigned j = 0; j < settings.numberStateVariableUpdates; j++) {
#pragma omp simd
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // fault strength using friction coefficient and fluid pressure from previous
        // timestep/iteration update state variable using sliprate from the previous time step
        localStateVariable[pointIndex] =
            static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                             ltsFace,
                                                             stateVarReference[pointIndex],
                                                             this->deltaT[timeIndex],
                                                             localSlipRate[pointIndex]);
      }
      this->tpMethod.calcFluidPressure(normalStress,
                                       this->mu,
                                       localSlipRate,
                                       this->deltaT[timeIndex],
                                       false,
                                       timeIndex,
                                       ltsFace);

      updateNormalStress(normalStress, faultStresses, timeIndex, ltsFace);

      // solve for new slip rate
      hasConverged = this->invertSlipRateIterative(
          ltsFace, localStateVariable, normalStress, absoluteShearStress, testSlipRate);

#pragma omp simd
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // update local slip rate, now using V=(Vnew+Vold)/2
        // For the next SV update, use the mean slip rate between the initial guess and the one
        // found (Kaneko 2008, step 6)
        localSlipRate[pointIndex] = 0.5 * (this->slipRateMagnitude[ltsFace][pointIndex] +
                                           std::fabs(testSlipRate[pointIndex]));

        // solve again for Vnew
        this->slipRateMagnitude[ltsFace][pointIndex] = std::fabs(testSlipRate[pointIndex]);

        // update friction coefficient based on new state variable and slip rate
        this->mu[ltsFace][pointIndex] =
            static_cast<Derived*>(this)->updateMu(ltsFace,
                                                  pointIndex,
                                                  this->slipRateMagnitude[ltsFace][pointIndex],
                                                  localStateVariable[pointIndex]);
      } // End of pointIndex-loop
    }
  }

  void calcSlipRateAndTraction(const std::array<real, misc::NumPaddedPoints>& stateVarReference,
                               const std::array<real, misc::NumPaddedPoints>& localSlipRate,
                               std::array<real, misc::NumPaddedPoints>& localStateVariable,
                               const std::array<real, misc::NumPaddedPoints>& normalStress,
                               const std::array<real, misc::NumPaddedPoints>& absoluteTraction,
                               const FaultStresses<Executor::Host>& faultStresses,
                               TractionResults<Executor::Host>& tractionResults,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // SV from mean slip rate in tmp
      localStateVariable[pointIndex] =
          static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                           ltsFace,
                                                           stateVarReference[pointIndex],
                                                           this->deltaT[timeIndex],
                                                           localSlipRate[pointIndex]);

      // update LocMu for next strength determination, only needed for last update
      this->mu[ltsFace][pointIndex] =
          static_cast<Derived*>(this)->updateMu(ltsFace,
                                                pointIndex,
                                                this->slipRateMagnitude[ltsFace][pointIndex],
                                                localStateVariable[pointIndex]);
      const real strength = -this->mu[ltsFace][pointIndex] * normalStress[pointIndex];
      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                                  faultStresses.traction1[timeIndex][pointIndex];
      const real totalTraction2 = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                                  faultStresses.traction2[timeIndex][pointIndex];

      const auto divisor =
          strength + this->impAndEta[ltsFace].etaS * this->slipRateMagnitude[ltsFace][pointIndex];
      this->slipRate1[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalTraction1 / divisor;
      this->slipRate2[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalTraction2 / divisor;

      // calculate traction
      tractionResults.traction1[timeIndex][pointIndex] =
          faultStresses.traction1[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * this->slipRate1[ltsFace][pointIndex];
      tractionResults.traction2[timeIndex][pointIndex] =
          faultStresses.traction2[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * this->slipRate2[ltsFace][pointIndex];
      this->traction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
      this->traction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];

      // Compute slip
      // ABS of locSlipRate removed as it would be the accumulated slip that is usually not needed
      // in the solver, see linear slip weakening
      this->accumulatedSlipMagnitude[ltsFace][pointIndex] +=
          this->slipRateMagnitude[ltsFace][pointIndex] * this->deltaT[timeIndex];

      // update directional slip
      this->slip1[ltsFace][pointIndex] +=
          this->slipRate1[ltsFace][pointIndex] * this->deltaT[timeIndex];
      this->slip2[ltsFace][pointIndex] +=
          this->slipRate2[ltsFace][pointIndex] * this->deltaT[timeIndex];
    }
  }

  void saveDynamicStressOutput(unsigned int face) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {

      if (this->ruptureTime[face][pointIndex] > 0.0 &&
          this->ruptureTime[face][pointIndex] <= this->mFullUpdateTime &&
          this->dynStressTimePending[face][pointIndex] &&
          this->mu[face][pointIndex] <=
              (this->drParameters->muW +
               0.05 * (this->drParameters->rsF0 - this->drParameters->muW))) {
        this->dynStressTime[face][pointIndex] = this->mFullUpdateTime;
        this->dynStressTimePending[face][pointIndex] = false;
      }
    }
  }

  /**
   * Solve for new slip rate (\f$\hat{s}\f$), applying the Newton-Raphson algorithm.
   * \f$\hat{s}\f$ has to fulfill
   * \f[g := \frac{1}{\eta_s} \cdot (\sigma_n \cdot \mu - \Theta) - \hat{s} = 0.\f] c.f. Carsten
   * Uphoff's dissertation eq. (4.57). Find root of \f$g\f$ with \f$g^\prime = \partial g / \partial
   * \hat{s}\f$: \f$\hat{s}_{i+1} = \hat{s}_i - ( g_i / g^\prime_i )\f$
   * @param ltsFace index of the face for which we invert the sliprate
   * @param localStateVariable \f$\psi\f$, needed to compute \f$\mu = f(\hat{s}, \psi)\f$
   * @param normalStress \f$\sigma_n\f$
   * @param absoluteShearStress \f$\Theta\f$
   * @param slipRateTest \f$\hat{s}\f$
   */
  bool invertSlipRateIterative(unsigned int ltsFace,
                               const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                               const std::array<real, misc::NumPaddedPoints>& normalStress,
                               const std::array<real, misc::NumPaddedPoints>& absoluteShearStress,
                               std::array<real, misc::NumPaddedPoints>& slipRateTest) {
    // Note that we need double precision here, since single precision led to NaNs.
    double muF[misc::NumPaddedPoints];
    double dMuF[misc::NumPaddedPoints];
    double g[misc::NumPaddedPoints];
    double dG[misc::NumPaddedPoints];

    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // first guess = sliprate value of the previous step
      slipRateTest[pointIndex] = this->slipRateMagnitude[ltsFace][pointIndex];
    }

    for (unsigned i = 0; i < settings.maxNumberSlipRateUpdates; i++) {
#pragma omp simd
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // calculate friction coefficient and objective function
        muF[pointIndex] = static_cast<Derived*>(this)->updateMu(
            ltsFace, pointIndex, slipRateTest[pointIndex], localStateVariable[pointIndex]);
        dMuF[pointIndex] = static_cast<Derived*>(this)->updateMuDerivative(
            ltsFace, pointIndex, slipRateTest[pointIndex], localStateVariable[pointIndex]);
        g[pointIndex] = -this->impAndEta[ltsFace].invEtaS *
                            (std::fabs(normalStress[pointIndex]) * muF[pointIndex] -
                             absoluteShearStress[pointIndex]) -
                        slipRateTest[pointIndex];
      }

      // max element of g must be smaller than newtonTolerance
      const bool hasConverged = std::all_of(std::begin(g), std::end(g), [&](auto val) {
        return std::fabs(val) < settings.newtonTolerance;
      });
      if (hasConverged) {
        return hasConverged;
      }
#pragma omp simd
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {

        // derivative of g
        dG[pointIndex] = -this->impAndEta[ltsFace].invEtaS *
                             (std::fabs(normalStress[pointIndex]) * dMuF[pointIndex]) -
                         1.0;
        // newton update
        const real tmp3 = g[pointIndex] / dG[pointIndex];
        slipRateTest[pointIndex] = std::max(rs::almostZero(), slipRateTest[pointIndex] - tmp3);
      }
    }
    return false;
  }

  void updateNormalStress(std::array<real, misc::NumPaddedPoints>& normalStress,
                          const FaultStresses<Executor::Host>& faultStresses,
                          size_t timeIndex,
                          size_t ltsFace) {
    // Todo(SW): consider poroelastic materials together with thermal pressurisation
#pragma omp simd
    for (size_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      normalStress[pointIndex] = std::min(static_cast<real>(0.0),
                                          faultStresses.normalStress[timeIndex][pointIndex] +
                                              this->initialStressInFaultCS[ltsFace][pointIndex][0] +
                                              faultStresses.fluidPressure[timeIndex][pointIndex] +
                                              this->initialPressure[ltsFace][pointIndex] -
                                              tpMethod.getFluidPressure(ltsFace, pointIndex));
    }
  }

  protected:
  // Attributes
  real (*a)[misc::NumPaddedPoints]{};
  real (*sl0)[misc::NumPaddedPoints]{};
  real (*stateVariable)[misc::NumPaddedPoints]{};

  TPMethod tpMethod;
  rs::Settings settings{};
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATE_H_
