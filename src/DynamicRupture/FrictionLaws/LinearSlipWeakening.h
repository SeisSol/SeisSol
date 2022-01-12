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

  /**
   * critical velocity at which slip rate is considered as being zero for instaneous healing
   */
  static constexpr real u_0 = 10e-14;

  real (*d_c)[numPaddedPoints];
  real (*mu_S)[numPaddedPoints];
  real (*mu_D)[numPaddedPoints];
  real (*cohesion)[numPaddedPoints];

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             std::array<real, numPaddedPoints>& stateVariableBuffer,
                             std::array<real, numPaddedPoints>& strengthBuffer,
                             unsigned int ltsFace,
                             unsigned int timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    this->calcStrengthHook(faultStresses, strengthBuffer, timeIndex, ltsFace);

    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    this->calcSlipRateAndTraction(faultStresses, strengthBuffer, timeIndex, ltsFace);

    // function g, output: stateVariablePsi & outputSlip
    this->calcStateVariableHook(stateVariableBuffer, timeIndex, ltsFace);

    // function f, output: calculated mu
    this->frictionFunctionHook(stateVariableBuffer, ltsFace);
  }

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
    this->d_c = layerData.var(concreteLts->d_c);
    this->mu_S = layerData.var(concreteLts->mu_s);
    this->mu_D = layerData.var(concreteLts->mu_d);
    this->cohesion = layerData.var(concreteLts->cohesion);
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slipStrike and slipDip
   */
  void calcSlipRateAndTraction(FaultStresses& faultStresses,
                               std::array<real, numPaddedPoints>& strength,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      real totalStressXY = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                           faultStresses.XYStressGP[timeIndex][pointIndex];
      real totalStressXZ = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                           faultStresses.XZStressGP[timeIndex][pointIndex];
      real absoluteShearStress = std::sqrt(std::pow(totalStressXY, 2) + std::pow(totalStressXZ, 2));

      // calculate slip rates
      this->slipRateMagnitude[ltsFace][pointIndex] = std::max(
          static_cast<real>(0.0),
          (absoluteShearStress - strength[pointIndex]) * this->impAndEta[ltsFace].inv_eta_s);

      this->slipRateStrike[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalStressXY / absoluteShearStress;
      this->slipRateDip[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalStressXZ / absoluteShearStress;

      // calculate traction
      faultStresses.XYTractionResultGP[timeIndex][pointIndex] =
          faultStresses.XYStressGP[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].eta_s * this->slipRateStrike[ltsFace][pointIndex];
      faultStresses.XZTractionResultGP[timeIndex][pointIndex] =
          faultStresses.XZStressGP[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].eta_s * this->slipRateDip[ltsFace][pointIndex];
      this->tractionXY[ltsFace][pointIndex] =
          faultStresses.XYTractionResultGP[timeIndex][pointIndex];
      this->tractionXZ[ltsFace][pointIndex] =
          faultStresses.XYTractionResultGP[timeIndex][pointIndex];

      // update directional slip
      this->slipStrike[ltsFace][pointIndex] +=
          this->slipRateStrike[ltsFace][pointIndex] * this->deltaT[timeIndex];
      this->slipDip[ltsFace][pointIndex] +=
          this->slipRateDip[ltsFace][pointIndex] * this->deltaT[timeIndex];
    }
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                            unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      this->mu[ltsFace][pointIndex] =
          mu_S[ltsFace][pointIndex] -
          (mu_S[ltsFace][pointIndex] - mu_D[ltsFace][pointIndex]) * stateVariablePsi[pointIndex];
    }
  }

  /*
   * instantaneous healing option Reset Mu and Slip
   */
  void instantaneousHealing(unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      if (this->slipRateMagnitude[ltsFace][pointIndex] < u_0) {
        this->mu[ltsFace][pointIndex] = mu_S[ltsFace][pointIndex];
        this->slip[ltsFace][pointIndex] = 0.0;
      }
    }
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  void saveDynamicStressOutput(unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

      if (this->ruptureTime[ltsFace][pointIndex] > 0.0 &&
          this->ruptureTime[ltsFace][pointIndex] <= this->m_fullUpdateTime &&
          this->ds[pointIndex] &&
          std::fabs(this->slip[ltsFace][pointIndex]) >= d_c[ltsFace][pointIndex]) {
        this->dynStressTime[ltsFace][pointIndex] = this->m_fullUpdateTime;
        this->ds[ltsFace][pointIndex] = false;
      }
    }
  }

  void calcStrengthHook(FaultStresses& faultStresses,
                        std::array<real, numPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace) {
    static_cast<Derived*>(this)->calcStrengthHook(faultStresses, strength, timeIndex, ltsFace);
  }

  void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace) {
    static_cast<Derived*>(this)->calcStateVariableHook(stateVariable, timeIndex, ltsFace);
  }
  void preHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace) {
    static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
  };
  void postHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace) {
    static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);
  };
};

class LinearSlipWeakeningLaw : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw> {
  public:
  using LinearSlipWeakeningBase::LinearSlipWeakeningBase;

  void preHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);
  void postHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);

  void calcStrengthHook(FaultStresses& faultStresses,
                        std::array<real, numPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace);

  void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace);
};

class LinearSlipWeakeningLawForcedRuptureTime : public LinearSlipWeakeningLaw {
  public:
  using LinearSlipWeakeningLaw::LinearSlipWeakeningLaw;

  real (*forcedRuptureTime)[numPaddedPoints];
  real* tn;

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  void preHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);
  void postHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace);

  void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariable,
                             unsigned int timeIndex,
                             unsigned int ltsFace);
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
                        std::array<real, numPaddedPoints>& strength,
                        unsigned int timeIndex,
                        unsigned int ltsFace);

  real (*regularisedStrength)[numPaddedPoints];

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  void prak_clif_mod(real& strength, real& sigma, real& LocSlipRate, real& mu, real& dt);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_LINEARSLIPWEAKENING_H
