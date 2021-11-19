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
  bool (*DS)[numPaddedPoints];
  real (*dynStressTime)[numPaddedPoints];

  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture* dynRup,
                real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                real fullUpdateTime,
                double timeWeights[CONVERGENCE_ORDER]) {

    BaseFrictionLaw<LinearSlipWeakeningBase<Derived>>::copyLtsTreeToLocal(
        layerData, dynRup, fullUpdateTime);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      // initialize struct for in/outputs stresses
      FaultStresses faultStresses = {};

      // declare local variables
      dynamicRupture::kernel::resampleParameter resampleKrnl;
      resampleKrnl.resampleM = init::resample::Values;

      std::array<real, numPaddedPoints> outputSlip{0};
      std::array<real, numPaddedPoints> stateVariablePsi{0};
      std::array<real, numPaddedPoints> strength{0};
      setTimeHook(ltsFace);

      this->precomputeStressFromQInterpolated(
          faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      for (int timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) { // loop over time steps
        // computes fault strength, which is the critical value whether active slip exists.
        static_cast<Derived*>(this)->calcStrengthHook(strength, faultStresses, timeIndex, ltsFace);

        // computes resulting slip rates, traction and slip dependent on current friction
        // coefficient and strength
        calcSlipRateAndTraction(strength, faultStresses, timeIndex, ltsFace);

        // function g, output: stateVariablePsi & outputSlip
        static_cast<Derived*>(this)->calcStateVariableHook(
            stateVariablePsi, outputSlip, resampleKrnl, timeIndex, ltsFace);

        // function f, output: calculated mu
        frictionFunctionHook(stateVariablePsi, ltsFace);

        // instantaneous healing option Reset Mu and Slip
        if (this->drParameters.isInstaHealingOn == true) {
          instantaneousHealing(ltsFace);
        }
      } // End of timeIndex-Loop

      // output rupture front
      this->saveRuptureFrontOutput(ltsFace);

      // output time when shear stress is equal to the dynamic stress after rupture arrived
      // currently only for linear slip weakening
      saveDynamicStressOutput(ltsFace);

      // output peak slip rate
      this->savePeakSlipRateOutput(ltsFace);

      //---compute and store slip to determine the magnitude of an earthquake ---
      //    to this end, here the slip is computed and averaged per element
      //    in calc_seissol.f90 this value will be multiplied by the element surface
      //    and an output happened once at the end of the simulation
      this->saveAverageSlipOutput(outputSlip, ltsFace);

      this->postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace],
                                                 QInterpolatedMinus[ltsFace],
                                                 faultStresses,
                                                 timeWeights,
                                                 ltsFace);
    } // End of Loop over Faces
  }   // End of Function evaluate

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
    this->d_c = layerData.var(concreteLts->d_c);
    this->mu_S = layerData.var(concreteLts->mu_s);
    this->mu_D = layerData.var(concreteLts->mu_d);
    this->cohesion = layerData.var(concreteLts->cohesion);
    this->DS = layerData.var(concreteLts->ds);
    this->averagedSlip = layerData.var(concreteLts->averagedSlip);
    this->dynStressTime = layerData.var(concreteLts->dynStressTime);
  }

  /**
   * Hook for FL16 to set tn equal to m_fullUpdateTime outside the calculation loop
   */
  void setTimeHook(unsigned int ltsFace) {}

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slipStrike and slipDip
   */
  void calcSlipRateAndTraction(std::array<real, numPaddedPoints>& strength,
                               FaultStresses& faultStresses,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
    std::array<real, numPaddedPoints> TotalShearStressYZ;
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      //-------------------------------------
      // calculate TotalShearStress in Y and Z direction
      TotalShearStressYZ[pointIndex] =
          std::sqrt(std::pow(this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                                 faultStresses.XYStressGP[timeIndex][pointIndex],
                             2) +
                    std::pow(this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                                 faultStresses.XZStressGP[timeIndex][pointIndex],
                             2));

      //-------------------------------------
      // calculate SlipRates
      this->slipRateMagnitude[ltsFace][pointIndex] =
          std::max(static_cast<real>(0.0),
                   (TotalShearStressYZ[pointIndex] - strength[pointIndex]) *
                       this->impAndEta[ltsFace].inv_eta_s);

      this->slipRateStrike[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] *
          (this->initialStressInFaultCS[ltsFace][pointIndex][3] +
           faultStresses.XYStressGP[timeIndex][pointIndex]) /
          TotalShearStressYZ[pointIndex];
      this->slipRateDip[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] *
          (this->initialStressInFaultCS[ltsFace][pointIndex][5] +
           faultStresses.XZStressGP[timeIndex][pointIndex]) /
          TotalShearStressYZ[pointIndex];

      //-------------------------------------
      // calculateTraction
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

      //-------------------------------------
      // update Directional Slip
      this->slipStrike[ltsFace][pointIndex] +=
          this->slipRateStrike[ltsFace][pointIndex] * this->deltaT[timeIndex];
      this->slipDip[ltsFace][pointIndex] +=
          this->slipRateDip[ltsFace][pointIndex] * this->deltaT[timeIndex];
    }
  }

  /*
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
          this->ruptureTime[ltsFace][pointIndex] <= this->m_fullUpdateTime && DS[pointIndex] &&
          std::fabs(this->slip[ltsFace][pointIndex]) >= d_c[ltsFace][pointIndex]) {
        dynStressTime[ltsFace][pointIndex] = this->m_fullUpdateTime;
        DS[ltsFace][pointIndex] = false;
      }
    }
  }
};

class LinearSlipWeakeningLaw : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw> {
  public:
  using LinearSlipWeakeningBase::LinearSlipWeakeningBase;
  void calcStrengthHook(std::array<real, numPaddedPoints>& strength,
                        FaultStresses& faultStresses,
                        unsigned int timeIndex,
                        unsigned int ltsFace);

  void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                             std::array<real, numPaddedPoints>& outputSlip,
                             dynamicRupture::kernel::resampleParameter& resampleKrnl,
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

  void setTimeHook(unsigned int ltsFace);

  void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                             std::array<real, numPaddedPoints>& outputSlip,
                             dynamicRupture::kernel::resampleParameter& resampleKrnl,
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
  void calcStrengthHook(std::array<real, numPaddedPoints>& strength,
                        FaultStresses& faultStresses,
                        unsigned int timeIndex,
                        unsigned int ltsFace);

  void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                             std::array<real, numPaddedPoints>& outputSlip,
                             dynamicRupture::kernel::resampleParameter& resampleKrnl,
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
