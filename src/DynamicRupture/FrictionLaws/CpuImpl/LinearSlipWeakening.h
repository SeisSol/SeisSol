// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_LINEARSLIPWEAKENING_H_

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
                             std::size_t ltsFace,
                             uint32_t timeIndex) {
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

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {
    this->dC = layerData.var<LTSLinearSlipWeakening::DC>();
    this->muS = layerData.var<LTSLinearSlipWeakening::MuS>();
    this->muD = layerData.var<LTSLinearSlipWeakening::MuD>();
    this->cohesion = layerData.var<LTSLinearSlipWeakening::Cohesion>();
    this->forcedRuptureTime = layerData.var<LTSLinearSlipWeakening::ForcedRuptureTime>();
    specialization.copyStorageToLocal(layerData);
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  void calcSlipRateAndTraction(const FaultStresses<Executor::Host>& faultStresses,
                               TractionResults<Executor::Host>& tractionResults,
                               std::array<real, misc::NumPaddedPoints>& strength,
                               uint32_t timeIndex,
                               std::size_t ltsFace) {
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 = this->initialStressInFaultCS[ltsFace][3][pointIndex] +
                                  faultStresses.traction1[timeIndex][pointIndex];
      const real totalTraction2 = this->initialStressInFaultCS[ltsFace][5][pointIndex] +
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

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {
  };
  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {
  };

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(std::array<real, misc::NumPaddedPoints>& stateVariable,
                            std::size_t ltsFace) {
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
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
  void saveDynamicStressOutput(std::size_t ltsFace, real time) {
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      if (this->dynStressTimePending[ltsFace][pointIndex] &&
          std::fabs(this->accumulatedSlipMagnitude[ltsFace][pointIndex]) >=
              dC[ltsFace][pointIndex]) {
        this->dynStressTime[ltsFace][pointIndex] = time;
        this->dynStressTimePending[ltsFace][pointIndex] = false;
      }
    }
  }

  void calcStrengthHook(const FaultStresses<Executor::Host>& faultStresses,
                        std::array<real, misc::NumPaddedPoints>& strength,
                        uint32_t timeIndex,
                        std::size_t ltsFace) {
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate fault strength (Uphoff eq 2.44) with addition cohesion term
      const real totalNormalStress = this->initialStressInFaultCS[ltsFace][0][pointIndex] +
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
                             uint32_t timeIndex,
                             std::size_t ltsFace) {
    alignas(Alignment) real resampledSlipRate[misc::NumPaddedPoints]{};
    specialization.resampleSlipRate(resampledSlipRate, this->slipRateMagnitude[ltsFace]);

    real time = this->mFullUpdateTime;
    for (uint32_t i = 0; i <= timeIndex; ++i) {
      time += this->deltaT[i];
    }
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
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
      if (this->drParameters->t0[0] == 0) {
        // avoid branching
        // if time > forcedRuptureTime, then f2 = 1.0, else f2 = 0.0
        f2 = 1.0 * static_cast<double>(time >= this->forcedRuptureTime[ltsFace][pointIndex]);
      } else {
        f2 = std::clamp((time - this->forcedRuptureTime[ltsFace][pointIndex]) /
                            this->drParameters->t0[0],
                        static_cast<real>(0.0),
                        static_cast<real>(1.0));
      }
      stateVariable[pointIndex] = std::max(stateVariable[pointIndex], f2);
    }
  }

  protected:
  real (*__restrict dC)[misc::NumPaddedPoints]{};
  real (*__restrict muS)[misc::NumPaddedPoints]{};
  real (*__restrict muD)[misc::NumPaddedPoints]{};
  real (*__restrict cohesion)[misc::NumPaddedPoints]{};
  real (*__restrict forcedRuptureTime)[misc::NumPaddedPoints]{};
  SpecializationT specialization;
};

class NoSpecialization {
  public:
  explicit NoSpecialization(seissol::initializer::parameters::DRParameters* parameters) {};

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {};
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
                                std::size_t ltsFace,
                                std::uint32_t pointIndex) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  }

#pragma omp declare simd
  static real strengthHook(real strength,
                           real localSlipRate,
                           real deltaT,
                           std::size_t ltsFace,
                           std::uint32_t pointIndex) {
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

  void copyStorageToLocal(DynamicRupture::Layer& layerData);
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
                                std::size_t ltsFace,
                                std::uint32_t pointIndex) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  }

#pragma omp declare simd
  real strengthHook(real strength,
                    real localSlipRate,
                    real deltaT,
                    std::size_t ltsFace,
                    std::uint32_t pointIndex);

  protected:
  seissol::initializer::parameters::DRParameters* drParameters;
  real (*__restrict regularizedStrength)[misc::NumPaddedPoints]{};
};

/**
 * Modified LSW friction as discussed in github issue #1058
 */
class TPApprox {
  public:
  explicit TPApprox(seissol::initializer::parameters::DRParameters* parameters)
      : drParameters(parameters) {};

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {}
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
                         std::size_t ltsFace,
                         std::uint32_t pointIndex);

#pragma omp declare simd
  static real strengthHook(real strength,
                           real localSlipRate,
                           real deltaT,
                           std::size_t ltsFace,
                           std::uint32_t pointIndex) {
    return strength;
  };

  protected:
  seissol::initializer::parameters::DRParameters* drParameters;
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_LINEARSLIPWEAKENING_H_
