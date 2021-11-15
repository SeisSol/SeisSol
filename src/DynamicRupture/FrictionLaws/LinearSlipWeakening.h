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
class LinearSlipWeakeningBase : public BaseFrictionLaw {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;

  protected:
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

  public:
  virtual void
      evaluate(seissol::initializers::Layer& layerData,
               seissol::initializers::DynamicRupture* dynRup,
               real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real fullUpdateTime,
               double timeWeights[CONVERGENCE_ORDER]) override {

    copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

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

      precomputeStressFromQInterpolated(
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
        if (drParameters.isInstaHealingOn == true) {
          instantaneousHealing(ltsFace);
        }
      } // End of timeIndex-Loop

      // output rupture front
      saveRuptureFrontOutput(ltsFace);

      // output time when shear stress is equal to the dynamic stress after rupture arrived
      // currently only for linear slip weakening
      saveDynamicStressOutput(ltsFace);

      // output peak slip rate
      savePeakSlipRateOutput(ltsFace);

      //---compute and store slip to determine the magnitude of an earthquake ---
      //    to this end, here the slip is computed and averaged per element
      //    in calc_seissol.f90 this value will be multiplied by the element surface
      //    and an output happened once at the end of the simulation
      saveAverageSlipOutput(outputSlip, ltsFace);

      postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace],
                                           QInterpolatedMinus[ltsFace],
                                           faultStresses,
                                           timeWeights,
                                           ltsFace);
    } // End of Loop over Faces
  }   // End of Function evaluate

  protected:
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) override {
    // first copy all Variables from the Base Lts dynRup tree
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
    d_c = layerData.var(concreteLts->d_c);
    mu_S = layerData.var(concreteLts->mu_s);
    mu_D = layerData.var(concreteLts->mu_d);
    cohesion = layerData.var(concreteLts->cohesion);
    DS = layerData.var(concreteLts->ds);
    averagedSlip = layerData.var(concreteLts->averagedSlip);
    dynStressTime = layerData.var(concreteLts->dynStressTime);
  }

  /**
   * Hook for FL16 to set tn equal to m_fullUpdateTime outside the calculation loop
   */
  virtual void setTimeHook(unsigned int ltsFace) {}

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slipStrike and slipDip
   */
  virtual void calcSlipRateAndTraction(std::array<real, numPaddedPoints>& strength,
                                       FaultStresses& faultStresses,
                                       unsigned int timeIndex,
                                       unsigned int ltsFace) {
    std::array<real, numPaddedPoints> TotalShearStressYZ;
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      //-------------------------------------
      // calculate TotalShearStress in Y and Z direction
      TotalShearStressYZ[pointIndex] =
          std::sqrt(std::pow(initialStressInFaultCS[ltsFace][pointIndex][3] +
                                 faultStresses.XYStressGP[timeIndex][pointIndex],
                             2) +
                    std::pow(initialStressInFaultCS[ltsFace][pointIndex][5] +
                                 faultStresses.XZStressGP[timeIndex][pointIndex],
                             2));

      //-------------------------------------
      // calculate SlipRates
      slipRateMagnitude[ltsFace][pointIndex] = std::max(
          static_cast<real>(0.0),
          (TotalShearStressYZ[pointIndex] - strength[pointIndex]) * impAndEta[ltsFace].inv_eta_s);

      slipRateStrike[ltsFace][pointIndex] = slipRateMagnitude[ltsFace][pointIndex] *
                                            (initialStressInFaultCS[ltsFace][pointIndex][3] +
                                             faultStresses.XYStressGP[timeIndex][pointIndex]) /
                                            TotalShearStressYZ[pointIndex];
      slipRateDip[ltsFace][pointIndex] = slipRateMagnitude[ltsFace][pointIndex] *
                                         (initialStressInFaultCS[ltsFace][pointIndex][5] +
                                          faultStresses.XZStressGP[timeIndex][pointIndex]) /
                                         TotalShearStressYZ[pointIndex];

      //-------------------------------------
      // calculateTraction
      faultStresses.XYTractionResultGP[timeIndex][pointIndex] =
          faultStresses.XYStressGP[timeIndex][pointIndex] -
          impAndEta[ltsFace].eta_s * slipRateStrike[ltsFace][pointIndex];
      faultStresses.XZTractionResultGP[timeIndex][pointIndex] =
          faultStresses.XZStressGP[timeIndex][pointIndex] -
          impAndEta[ltsFace].eta_s * slipRateDip[ltsFace][pointIndex];
      tractionXY[ltsFace][pointIndex] = faultStresses.XYTractionResultGP[timeIndex][pointIndex];
      tractionXZ[ltsFace][pointIndex] = faultStresses.XYTractionResultGP[timeIndex][pointIndex];

      //-------------------------------------
      // update Directional Slip
      slipStrike[ltsFace][pointIndex] += slipRateStrike[ltsFace][pointIndex] * deltaT[timeIndex];
      slipDip[ltsFace][pointIndex] += slipRateDip[ltsFace][pointIndex] * deltaT[timeIndex];
    }
  }

  /*
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  virtual void frictionFunctionHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                                    unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      mu[ltsFace][pointIndex] =
          mu_S[ltsFace][pointIndex] -
          (mu_S[ltsFace][pointIndex] - mu_D[ltsFace][pointIndex]) * stateVariablePsi[pointIndex];
    }
  }

  /*
   * instantaneous healing option Reset Mu and Slip
   */
  virtual void instantaneousHealing(unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      if (slipRateMagnitude[ltsFace][pointIndex] < u_0) {
        mu[ltsFace][pointIndex] = mu_S[ltsFace][pointIndex];
        slip[ltsFace][pointIndex] = 0.0;
      }
    }
  }
  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  virtual void saveDynamicStressOutput(unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

      if (ruptureTime[ltsFace][pointIndex] > 0.0 &&
          ruptureTime[ltsFace][pointIndex] <= m_fullUpdateTime && DS[pointIndex] &&
          std::fabs(slip[ltsFace][pointIndex]) >= d_c[ltsFace][pointIndex]) {
        dynStressTime[ltsFace][pointIndex] = m_fullUpdateTime;
        DS[ltsFace][pointIndex] = false;
      }
    }
  }
};

class LinearSlipWeakeningLaw : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw> {
  public:
  using LinearSlipWeakeningBase::LinearSlipWeakeningBase;
  virtual void calcStrengthHook(std::array<real, numPaddedPoints>& strength,
                                FaultStresses& faultStresses,
                                unsigned int timeIndex,
                                unsigned int ltsFace);

  virtual void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                                     std::array<real, numPaddedPoints>& outputSlip,
                                     dynamicRupture::kernel::resampleParameter& resampleKrnl,
                                     unsigned int timeIndex,
                                     unsigned int ltsFace);
};

class LinearSlipWeakeningLawForcedRuptureTime : public LinearSlipWeakeningLaw {
  public:
  using LinearSlipWeakeningLaw::LinearSlipWeakeningLaw;

  protected:
  real (*forcedRuptureTime)[numPaddedPoints];
  real* tn;

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) override;

  virtual void setTimeHook(unsigned int ltsFace) override;

  virtual void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                                     std::array<real, numPaddedPoints>& outputSlip,
                                     dynamicRupture::kernel::resampleParameter& resampleKrnl,
                                     unsigned int timeIndex,
                                     unsigned int ltsFace) override;
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
  virtual void calcStrengthHook(std::array<real, numPaddedPoints>& strength,
                                FaultStresses& faultStresses,
                                unsigned int timeIndex,
                                unsigned int ltsFace);

  virtual void calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariablePsi,
                                     std::array<real, numPaddedPoints>& outputSlip,
                                     dynamicRupture::kernel::resampleParameter& resampleKrnl,
                                     unsigned int timeIndex,
                                     unsigned int ltsFace);

  protected:
  real (*regularisedStrength)[numPaddedPoints];

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) override;

  void prak_clif_mod(real& strength, real& sigma, real& LocSlipRate, real& mu, real& dt);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_LINEARSLIPWEAKENING_H
