#ifndef SEISSOL_LINEARSLIPWEAKENING_H
#define SEISSOL_LINEARSLIPWEAKENING_H

#include "BaseFrictionLaw.h"

#include "utils/logger.h"

namespace seissol::dr::friction_law {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <typename Config, class SpecializationT>
class LinearSlipWeakeningLaw
    : public BaseFrictionLaw<Config, LinearSlipWeakeningLaw<Config, SpecializationT>> {
  public:
  using RealT = typename Config::RealT;
  explicit LinearSlipWeakeningLaw(DRParameters* drParameters)
      : BaseFrictionLaw<Config, LinearSlipWeakeningLaw<Config, SpecializationT>>(drParameters),
        specialization(drParameters) {}

  void updateFrictionAndSlip(FaultStresses<Config> const& faultStresses,
                             TractionResults<Config>& tractionResults,
                             std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
                             std::array<RealT, misc::numPaddedPoints<Config>>& strengthBuffer,
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

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime) {
    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSLinearSlipWeakening<Config> const* const>(dynRup);
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
  void calcSlipRateAndTraction(FaultStresses<Config> const& faultStresses,
                               TractionResults<Config>& tractionResults,
                               std::array<RealT, misc::numPaddedPoints<Config>>& strength,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      const RealT totalTraction1 = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                                   faultStresses.traction1[timeIndex][pointIndex];
      const RealT totalTraction2 = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                                   faultStresses.traction2[timeIndex][pointIndex];
      const RealT absoluteTraction = misc::magnitude(totalTraction1, totalTraction2);

      // calculate slip rates
      this->slipRateMagnitude[ltsFace][pointIndex] =
          std::max(static_cast<RealT>(0.0),
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

  void preHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
               unsigned int ltsFace){};
  void postHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariableBuffer,
                unsigned int ltsFace){};

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariable,
                            unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
      this->mu[ltsFace][pointIndex] =
          muS[ltsFace][pointIndex] -
          (muS[ltsFace][pointIndex] - muD[ltsFace][pointIndex]) * stateVariable[pointIndex];
    }
  }

  /**
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   */
  void saveDynamicStressOutput(unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
      if (this->dynStressTimePending[ltsFace][pointIndex] &&
          std::fabs(this->accumulatedSlipMagnitude[ltsFace][pointIndex]) >=
              dC[ltsFace][pointIndex]) {
        this->dynStressTime[ltsFace][pointIndex] = this->mFullUpdateTime;
        this->dynStressTimePending[ltsFace][pointIndex] = false;
      }
    }
  }

  void calcStrengthHook(FaultStresses<Config> const& faultStresses,
                        std::array<RealT, misc::numPaddedPoints<Config>>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace) {
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
      // calculate fault strength (Uphoff eq 2.44) with addition cohesion term
      const RealT totalNormalStress = this->initialStressInFaultCS[ltsFace][pointIndex][0] +
                                      faultStresses.normalStress[timeIndex][pointIndex];
      strength[pointIndex] =
          -cohesion[ltsFace][pointIndex] -
          this->mu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<RealT>(0.0));

      strength[pointIndex] =
          specialization.strengthHook(strength[pointIndex],
                                      this->slipRateMagnitude[ltsFace][pointIndex],
                                      this->deltaT[timeIndex],
                                      ltsFace,
                                      pointIndex);
    }
  }

  void calcStateVariableHook(std::array<RealT, misc::numPaddedPoints<Config>>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace) {
    alignas(Alignment) RealT resampledSlipRate[misc::numPaddedPoints<Config>]{};
    specialization.resampleSlipRate(resampledSlipRate, this->slipRateMagnitude[ltsFace]);

    const RealT time = this->mFullUpdateTime + this->deltaT[timeIndex];
#pragma omp simd
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
      // integrate slip rate to get slip = state variable

      const auto update = resampledSlipRate[pointIndex] * this->deltaT[timeIndex];
      this->accumulatedSlipMagnitude[ltsFace][pointIndex] += update;

      // Actually slip is already the stateVariable for this FL, but to simplify the next equations
      // we divide it here by the critical distance.
      stateVariable[pointIndex] = std::min(
          std::fabs(this->accumulatedSlipMagnitude[ltsFace][pointIndex]) / dC[ltsFace][pointIndex],
          static_cast<RealT>(1.0));

      // Forced rupture time
      RealT f2 = 0.0;
      if (this->drParameters->t0 == 0) {
        // avoid branching
        // if time > forcedRuptureTime, then f2 = 1.0, else f2 = 0.0
        f2 = 1.0 * (time >= this->forcedRuptureTime[ltsFace][pointIndex]);
      } else {
        f2 = std::clamp((time - this->forcedRuptureTime[ltsFace][pointIndex]) /
                            this->drParameters->t0,
                        static_cast<RealT>(0.0),
                        static_cast<RealT>(1.0));
      }
      stateVariable[pointIndex] = std::max(stateVariable[pointIndex], f2);
    }
  }

  protected:
  RealT (*dC)[misc::numPaddedPoints<Config>];
  RealT (*muS)[misc::numPaddedPoints<Config>];
  RealT (*muD)[misc::numPaddedPoints<Config>];
  RealT (*cohesion)[misc::numPaddedPoints<Config>];
  RealT (*forcedRuptureTime)[misc::numPaddedPoints<Config>];
  SpecializationT specialization;
};

template <typename Config>
class NoSpecialization {
  public:
  using RealT = typename Config::RealT;
  explicit NoSpecialization(DRParameters* parameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime){};
  /**
   * Resample slip-rate, such that the state increment (slip) lies in the same polynomial space as
   * the degrees of freedom resampleMatrix first projects LocSR on the two-dimensional basis on
   * the reference triangle with degree less or equal than ConvergenceOrder-1, and then evaluates
   * the polynomial at the quadrature points
   */
  void resampleSlipRate(RealT (&resampledSlipRate)[dr::misc::numPaddedPoints<Config>],
                        RealT const (&slipRate)[dr::misc::numPaddedPoints<Config>]);
#pragma omp declare simd
  RealT strengthHook(RealT strength,
                     RealT localSlipRate,
                     RealT deltaT,
                     unsigned int ltsFace,
                     unsigned int pointIndex) {
    return strength;
  };
};

/**
 * Law for bimaterial faults, implements strength regularization (according to Prakash-Clifton)
 */
template <typename Config>
class BiMaterialFault {
  public:
  using RealT = typename Config::RealT;
  explicit BiMaterialFault(DRParameters* parameters) : drParameters(parameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime);
  /**
   * Resampling of the sliprate introduces artificial oscillations into the solution, if we use it
   * together with Prakash-Clifton regularization, so for the BiMaterialFault specialization, we
   * replace the resampling with a simple copy.
   */
  void resampleSlipRate(RealT (&resampledSlipRate)[dr::misc::numPaddedPoints<Config>],
                        RealT const (&slipRate)[dr::misc::numPaddedPoints<Config>]) {
    std::copy(std::begin(slipRate), std::end(slipRate), std::begin(resampledSlipRate));
  };
#pragma omp declare simd
  RealT strengthHook(RealT strength,
                     RealT localSlipRate,
                     RealT deltaT,
                     unsigned int ltsFace,
                     unsigned int pointIndex);

  protected:
  DRParameters* drParameters;
  RealT (*regularisedStrength)[misc::numPaddedPoints<Config>];
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_LINEARSLIPWEAKENING_H
