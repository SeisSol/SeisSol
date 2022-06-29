#ifndef SEISSOL_LINEARSLIPWEAKENING_H
#define SEISSOL_LINEARSLIPWEAKENING_H

#include "BaseFrictionLaw.h"

#include "utils/logger.h"

namespace seissol::dr::friction_law {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <class Derived>
class LinearSlipWeakeningBase : public BaseFrictionLaw<LinearSlipWeakeningBase<Derived>> {
  public:
  using BaseFrictionLaw<LinearSlipWeakeningBase<Derived>>::BaseFrictionLaw;

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             TractionResults& tractionResults,
                             std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::numPaddedPoints>& strengthBuffer,
                             unsigned int ltsFace,
                             unsigned int timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    this->calcStrengthHook(faultStresses, strengthBuffer, timeIndex, ltsFace);

    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    this->calcSlipRateAndTraction(
        faultStresses, tractionResults, strengthBuffer, timeIndex, ltsFace);

    this->calcStateVariableHook(stateVariableBuffer, timeIndex, ltsFace);

    this->frictionFunctionHook(stateVariableBuffer, ltsFace);
  }

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
    this->dC = layerData.var(concreteLts->dC);
    this->muS = layerData.var(concreteLts->muS);
    this->muD = layerData.var(concreteLts->muD);
    this->cohesion = layerData.var(concreteLts->cohesion);
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  void calcSlipRateAndTraction(FaultStresses& faultStresses,
                               TractionResults& tractionResults,
                               std::array<real, misc::numPaddedPoints>& strength,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      real totalTraction1 = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                            faultStresses.traction1[timeIndex][pointIndex];
      real totalTraction2 = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                            faultStresses.traction2[timeIndex][pointIndex];
      real absoluteTraction = misc::magnitude(totalTraction1, totalTraction2);

      // calculate slip rates
      this->slipRateMagnitude[ltsFace][pointIndex] =
          std::max(static_cast<real>(0.0),
                   (absoluteTraction - strength[pointIndex]) * this->impAndEta[ltsFace].invEtaS);

      auto divisor = strength[pointIndex] +
                     this->impAndEta[ltsFace].etaS * this->slipRateMagnitude[ltsFace][pointIndex];
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

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                            unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      this->mu[ltsFace][pointIndex] =
          muS[ltsFace][pointIndex] -
          (muS[ltsFace][pointIndex] - muD[ltsFace][pointIndex]) * stateVariable[pointIndex];
    }
  }

  /**
   * Instantaneous healing option:
   * Reset Mu and Slip, if slipRateMagnitude drops below threshold
   * This function is currently not used, as we miss an appropriate benchmark.
   */
  void instantaneousHealing(unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      if (this->slipRateMagnitude[ltsFace][pointIndex] < u0) {
        this->mu[ltsFace][pointIndex] = muS[ltsFace][pointIndex];
        this->accumulatedSlipMagnitude[ltsFace][pointIndex] = 0.0;
      }
    }
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  void saveDynamicStressOutput(unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {

      if (this->dynStressTimePending[ltsFace][pointIndex] &&
          std::fabs(this->accumulatedSlipMagnitude[ltsFace][pointIndex]) >=
              dC[ltsFace][pointIndex]) {
        this->dynStressTime[ltsFace][pointIndex] = this->mFullUpdateTime;
        this->dynStressTimePending[ltsFace][pointIndex] = false;
      }
    }
  }

  void calcStrengthHook(FaultStresses& faultStresses,
                        std::array<real, misc::numPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace) {
    static_cast<Derived*>(this)->calcStrengthHook(faultStresses, strength, timeIndex, ltsFace);
  }

  void calcStateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace) {
    static_cast<Derived*>(this)->calcStateVariableHook(stateVariable, timeIndex, ltsFace);
  }
  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace) {
    static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
  }
  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                unsigned int ltsFace) {
    static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);
  }

  protected:
  /**
   * critical velocity at which slip rate is considered as being zero for instaneous healing
   */
  static constexpr real u0 = 10e-14;
  real (*dC)[misc::numPaddedPoints];
  real (*muS)[misc::numPaddedPoints];
  real (*muD)[misc::numPaddedPoints];
  real (*cohesion)[misc::numPaddedPoints];
};

class LinearSlipWeakeningLaw : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw> {
  public:
  using LinearSlipWeakeningBase::LinearSlipWeakeningBase;

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);
  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);

  void calcStrengthHook(FaultStresses& faultStresses,
                        std::array<real, misc::numPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace);

  void calcStateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace);
};

class LinearSlipWeakeningLawForcedRuptureTime : public LinearSlipWeakeningLaw {
  public:
  using LinearSlipWeakeningLaw::LinearSlipWeakeningLaw;

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);
  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);

  void calcStateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace);

  protected:
  real (*forcedRuptureTime)[misc::numPaddedPoints];
  real* tn;
};

/**
 * Law for Bimaterial faults, implements strength regularization (according to prakash clifton)
 * currently regularized strength is not used (bug)
 * State variable (slip) is not resampled in this friction law!
 */
class LinearSlipWeakeningLawBimaterial
    : public LinearSlipWeakeningBase<LinearSlipWeakeningLawBimaterial> {
  public:
  using LinearSlipWeakeningBase::LinearSlipWeakeningBase;
  void calcStrengthHook(FaultStresses& faultStresses,
                        std::array<real, misc::numPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace);

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
               unsigned int ltsFace){};
  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                unsigned int ltsFace){};

  void calcStateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  protected:
  real (*regularisedStrength)[misc::numPaddedPoints];
  void prak_clif_mod(real& strength, real& sigma, real& locSlipRate, real& mu, real& dt);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_LINEARSLIPWEAKENING_H
