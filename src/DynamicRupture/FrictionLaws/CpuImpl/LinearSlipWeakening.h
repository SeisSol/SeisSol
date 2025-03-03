// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_LINEARSLIPWEAKENING_H_

#include "BaseFrictionLaw.h"

#include "utils/logger.h"

namespace seissol::dr::friction_law::cpu {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <class SpecializationT>
class LinearSlipWeakeningLaw : public BaseFrictionLaw<LinearSlipWeakeningLaw<SpecializationT>> {
  public:
  explicit LinearSlipWeakeningLaw(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionLaw<LinearSlipWeakeningLaw<SpecializationT>>(drParameters),
        specialization(drParameters) {}

  void updateFrictionAndSlip(const FaultStresses<Executor::Host>& faultStresses,
                             TractionResults<Executor::Host>& tractionResults,
                             std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::NumPaddedPoints>& strengthBuffer,
                             unsigned int ltsFace,
                             unsigned int timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    this->calcStrengthHook(faultStresses, strengthBuffer, timeIndex, ltsFace);

    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    this->calcSlipRateAndTraction(
        faultStresses, tractionResults, strengthBuffer, timeIndex, ltsFace);

    // integrate state variable in time
    this->calcStateVariableHook(stateVariableBuffer, timeIndex, ltsFace);

    // compute friction coefficient based on state variable and slip
    this->frictionFunctionHook(stateVariableBuffer, ltsFace);
  }

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening*>(dynRup);
    this->dC = layerData.var(concreteLts->dC);
    this->muS = layerData.var(concreteLts->muS);
    this->muD = layerData.var(concreteLts->muD);
    this->cohesion = layerData.var(concreteLts->cohesion);
    this->forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime);
    specialization.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  void calcSlipRateAndTraction(const FaultStresses<Executor::Host>& faultStresses,
                               TractionResults<Executor::Host>& tractionResults,
                               std::array<real, misc::NumPaddedPoints>& strength,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                                  faultStresses.traction1[timeIndex][pointIndex];
      const real totalTraction2 = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                                  faultStresses.traction2[timeIndex][pointIndex];
      const real absoluteTraction = misc::magnitude(totalTraction1, totalTraction2);

      // calculate slip rates
      this->slipRateMagnitude[ltsFace][pointIndex] =
          std::max(static_cast<real>(0.0),
                   (absoluteTraction - strength[pointIndex]) * this->impAndEta[ltsFace].invEtaS);

      const auto divisor = strength[pointIndex] + this->impAndEta[ltsFace].etaS *
                                                      this->slipRateMagnitude[ltsFace][pointIndex];
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

      // update directional slip
      this->slip1[ltsFace][pointIndex] +=
          this->slipRate1[ltsFace][pointIndex] * this->deltaT[timeIndex];
      this->slip2[ltsFace][pointIndex] +=
          this->slipRate2[ltsFace][pointIndex] * this->deltaT[timeIndex];
    }
  }

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, unsigned int ltsFace) {
  };
  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                unsigned int ltsFace) {};

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(std::array<real, misc::NumPaddedPoints>& stateVariable,
                            unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->mu[ltsFace][pointIndex] =
          muS[ltsFace][pointIndex] -
          (muS[ltsFace][pointIndex] - muD[ltsFace][pointIndex]) * stateVariable[pointIndex];
      // instantaneous healing
      if ((this->peakSlipRate[ltsFace][pointIndex] > this->drParameters->healingThreshold) &&
          (this->slipRateMagnitude[ltsFace][pointIndex] < this->drParameters->healingThreshold)) {
        this->mu[ltsFace][pointIndex] = muS[ltsFace][pointIndex];
        stateVariable[pointIndex] = 0.0;
      }
    }
  }

  /**
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   */
  void saveDynamicStressOutput(unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      if (this->dynStressTimePending[ltsFace][pointIndex] &&
          std::fabs(this->accumulatedSlipMagnitude[ltsFace][pointIndex]) >=
              dC[ltsFace][pointIndex]) {
        this->dynStressTime[ltsFace][pointIndex] = this->mFullUpdateTime;
        this->dynStressTimePending[ltsFace][pointIndex] = false;
      }
    }
  }

  void calcStrengthHook(const FaultStresses<Executor::Host>& faultStresses,
                        std::array<real, misc::NumPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate fault strength (Uphoff eq 2.44) with addition cohesion term
      const real totalNormalStress = this->initialStressInFaultCS[ltsFace][pointIndex][0] +
                                     faultStresses.normalStress[timeIndex][pointIndex] +
                                     this->initialPressure[ltsFace][pointIndex] +
                                     faultStresses.fluidPressure[timeIndex][pointIndex];

      strength[pointIndex] =
          -cohesion[ltsFace][pointIndex] -
          this->mu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<real>(0.0));

      strength[pointIndex] =
          specialization.strengthHook(strength[pointIndex],
                                      this->slipRateMagnitude[ltsFace][pointIndex],
                                      this->deltaT[timeIndex],
                                      ltsFace,
                                      pointIndex);
    }
  }

  void calcStateVariableHook(std::array<real, misc::NumPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace) {
    alignas(Alignment) real resampledSlipRate[misc::NumPaddedPoints]{};
    specialization.resampleSlipRate(resampledSlipRate, this->slipRateMagnitude[ltsFace]);

    const real time = this->mFullUpdateTime + this->deltaT[timeIndex];
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // integrate slip rate to get slip = state variable

      const auto update = resampledSlipRate[pointIndex] * this->deltaT[timeIndex];
      this->accumulatedSlipMagnitude[ltsFace][pointIndex] += update;

      // Actually slip is already the stateVariable for this FL, but to simplify the next equations
      // we divide it here by the critical distance.
      stateVariable[pointIndex] =
          specialization.stateVariableHook(this->accumulatedSlipMagnitude[ltsFace][pointIndex],
                                           dC[ltsFace][pointIndex],
                                           ltsFace,
                                           pointIndex);

      // Forced rupture time
      real f2 = 0.0;
      if (this->drParameters->t0 == 0) {
        // avoid branching
        // if time > forcedRuptureTime, then f2 = 1.0, else f2 = 0.0
        f2 = 1.0 * static_cast<double>(time >= this->forcedRuptureTime[ltsFace][pointIndex]);
      } else {
        f2 = std::clamp((time - this->forcedRuptureTime[ltsFace][pointIndex]) /
                            this->drParameters->t0,
                        static_cast<real>(0.0),
                        static_cast<real>(1.0));
      }
      stateVariable[pointIndex] = std::max(stateVariable[pointIndex], f2);
    }
  }

  protected:
  real (*dC)[misc::NumPaddedPoints]{};
  real (*muS)[misc::NumPaddedPoints]{};
  real (*muD)[misc::NumPaddedPoints]{};
  real (*cohesion)[misc::NumPaddedPoints]{};
  real (*forcedRuptureTime)[misc::NumPaddedPoints]{};
  SpecializationT specialization;
};

class NoSpecialization {
  public:
  explicit NoSpecialization(seissol::initializer::parameters::DRParameters* parameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {};
  /**
   * Resample slip-rate, such that the state increment (slip) lies in the same polynomial space as
   * the degrees of freedom resampleMatrix first projects LocSR on the two-dimensional basis on
   * the reference triangle with degree less or equal than ConvergenceOrder-1, and then evaluates
   * the polynomial at the quadrature points
   */
  static void resampleSlipRate(real (&resampledSlipRate)[dr::misc::NumPaddedPoints],
                               const real (&slipRate)[dr::misc::NumPaddedPoints]);
#pragma omp declare simd
  static real stateVariableHook(real localAccumulatedSlip,
                                real localDc,
                                unsigned int ltsFace,
                                unsigned int pointIndex) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  }

#pragma omp declare simd
  static real strengthHook(real strength,
                           real localSlipRate,
                           real deltaT,
                           unsigned int ltsFace,
                           unsigned int pointIndex) {
    return strength;
  };
};

/**
 * Law for bimaterial faults, implements strength regularization (according to Prakash-Clifton)
 */
class BiMaterialFault {
  public:
  explicit BiMaterialFault(seissol::initializer::parameters::DRParameters* parameters)
      : drParameters(parameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* dynRup,
                          real fullUpdateTime);
  /**
   * Resampling of the sliprate introduces artificial oscillations into the solution, if we use it
   * together with Prakash-Clifton regularization, so for the BiMaterialFault specialization, we
   * replace the resampling with a simple copy.
   */
  static void resampleSlipRate(real (&resampledSlipRate)[dr::misc::NumPaddedPoints],
                               const real (&slipRate)[dr::misc::NumPaddedPoints]) {
    std::copy(std::begin(slipRate), std::end(slipRate), std::begin(resampledSlipRate));
  };

#pragma omp declare simd
  static real stateVariableHook(real localAccumulatedSlip,
                                real localDc,
                                unsigned int ltsFace,
                                unsigned int pointIndex) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  }

#pragma omp declare simd
  real strengthHook(real strength,
                    real localSlipRate,
                    real deltaT,
                    unsigned int ltsFace,
                    unsigned int pointIndex);

  protected:
  seissol::initializer::parameters::DRParameters* drParameters;
  real (*regularizedStrength)[misc::NumPaddedPoints]{};
};

/**
 * Modified LSW friction as discussed in github issue #1058
 */
class TPApprox {
  public:
  explicit TPApprox(seissol::initializer::parameters::DRParameters* parameters)
      : drParameters(parameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}
  /**
   * Use a simple copy for now, maybe use proper resampling later
   */
  static void resampleSlipRate(real (&resampledSlipRate)[dr::misc::NumPaddedPoints],
                               const real (&slipRate)[dr::misc::NumPaddedPoints]) {
    std::copy(std::begin(slipRate), std::end(slipRate), std::begin(resampledSlipRate));
  };

#pragma omp declare simd
  real stateVariableHook(real localAccumulatedSlip,
                         real localDc,
                         unsigned int ltsFace,
                         unsigned int pointIndex);

#pragma omp declare simd
  static real strengthHook(real strength,
                           real localSlipRate,
                           real deltaT,
                           unsigned int ltsFace,
                           unsigned int pointIndex) {
    return strength;
  };

  protected:
  seissol::initializer::parameters::DRParameters* drParameters;
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_LINEARSLIPWEAKENING_H_
