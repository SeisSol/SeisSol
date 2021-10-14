#ifndef SEISSOL_LINEARSLIPWEAKENING_H
#define SEISSOL_LINEARSLIPWEAKENING_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
  template <class Derived>
  class LinearSlipWeakeningLaw;      //generalization of Linear slip weakening laws
  class LinearSlipWeakeningLawFL2;   //linear slip weakening
  class LinearSlipWeakeningLawFL16;  //Linear slip weakening forced time rapture
  class LinearSlipWeakeningLawBimaterialFL6;  //solver for bimaterial faults, currently has a bug, solution to the bug inside the function
};

/*
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
  template <class Derived>
  class seissol::dr::friction_law::LinearSlipWeakeningLaw : public seissol::dr::friction_law::BaseFrictionLaw {
  protected:
    static constexpr real   u_0  = 10e-14;   //critical velocity at which slip rate is considered as being zero for instaneous healing
    real                    (*d_c)[numOfPointsPadded];
    real                    (*mu_S)[numOfPointsPadded];
    real                    (*mu_D)[numOfPointsPadded];
    bool                    (*DS)[numOfPointsPadded];
    real                    (*dynStress_time)[numOfPointsPadded];

  public:
    virtual void evaluate(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup,
                          real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real fullUpdateTime,
                          double timeWeights[CONVERGENCE_ORDER]) override {

      copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
        //initialize struct for in/outputs stresses
        FaultStresses faultStresses = {};

        //declare local variables
        dynamicRupture::kernel::resampleParameter resampleKrnl;
        resampleKrnl.resampleM = init::resample::Values;

        std::array<real, numOfPointsPadded> outputSlip{0};
        std::array<real, numOfPointsPadded> stateVariablePsi{0};
        std::array<real, numOfPointsPadded> Strength{0};
        setTimeHook(ltsFace);

        precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

        for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
          //computes fault strength, which is the critical value whether active slip exists.
          static_cast<Derived*>(this)->calcStrengthHook(Strength, faultStresses, iTimeGP, ltsFace);

          //computes resulting slip rates, traction and slip dependent on current friction coefficient and Strength
          calcSlipRateAndTraction(Strength, faultStresses, iTimeGP , ltsFace);

          //function g, output: stateVariablePsi & outputSlip
          static_cast<Derived*>(this)->calcStateVariableHook(stateVariablePsi, outputSlip, resampleKrnl, iTimeGP, ltsFace);

          //function f, output: calculated mu
          frictionFunctionHook(stateVariablePsi, ltsFace);

          //instantaneous healing option Reset Mu and Slip
          if (m_Params->IsInstaHealingOn == true) {
            instantaneousHealing(ltsFace);
          }
        }//End of iTimeGP-Loop

        // output rupture front
        saveRuptureFrontOutput(ltsFace);

        //output time when shear stress is equal to the dynamic stress after rupture arrived
        //currently only for linear slip weakening
        saveDynamicStressOutput(ltsFace);

        //output peak slip rate
        savePeakSlipRateOutput(ltsFace);

        //---compute and store slip to determine the magnitude of an earthquake ---
        //    to this end, here the slip is computed and averaged per element
        //    in calc_seissol.f90 this value will be multiplied by the element surface
        //    and an output happened once at the end of the simulation
        saveAverageSlipOutput(outputSlip, ltsFace);

        postcomputeImposedStateFromNewStress(
            QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace],
            faultStresses, timeWeights, ltsFace);
      }//End of Loop over Faces
    }//End of Function evaluate

  protected:
    /*
     * copies all parameters from the DynamicRupture LTS to the local attributes
     */
    void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                            seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) override {
      //first copy all Variables from the Base Lts dynRup tree
      BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

      seissol::initializers::LTS_LinearSlipWeakeningFL2 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningFL2 *>(dynRup);
      d_c                                           = layerData.var(ConcreteLts->d_c);
      mu_S                                          = layerData.var(ConcreteLts->mu_S);
      mu_D                                          = layerData.var(ConcreteLts->mu_D);
      DS                                            = layerData.var(ConcreteLts->DS);
      averaged_Slip                                 = layerData.var(ConcreteLts->averaged_Slip);
      dynStress_time                                = layerData.var(ConcreteLts->dynStress_time);
    }


    /*
     * Hook for FL16 to set tn equal to m_fullUpdateTime outside the calculation loop
     */
    virtual void setTimeHook(unsigned int ltsFace) {}

    /*
     *  compute the slip rate and the traction from the fault strength and fault stresses
     *  also updates the directional slipStrike and slipDip
     */
    virtual void calcSlipRateAndTraction(
        std::array<real, numOfPointsPadded> &Strength,
        FaultStresses &faultStresses,
        unsigned int iTimeGP,
        unsigned int ltsFace
    ){
      std::array<real, numOfPointsPadded> TotalShearStressYZ;
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
        //-------------------------------------
        //calculate TotalShearStress in Y and Z direction
        TotalShearStressYZ[iBndGP] = std::sqrt(
            seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP], 2) +
            seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP], 2));

        //-------------------------------------
        // calculate SlipRates
        SlipRateMagnitude[ltsFace][iBndGP] = std::max((real)0.0, (TotalShearStressYZ[iBndGP] - Strength[iBndGP]) * impAndEta[ltsFace].inv_eta_s);

        slipRateStrike[ltsFace][iBndGP] = SlipRateMagnitude[ltsFace][iBndGP] * (initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP]) / TotalShearStressYZ[iBndGP];
        slipRateDip[ltsFace][iBndGP]  = SlipRateMagnitude[ltsFace][iBndGP] * (initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP]) / TotalShearStressYZ[iBndGP];

        //-------------------------------------
        //calculateTraction
        faultStresses.XYTractionResultGP[iTimeGP][iBndGP] = faultStresses.XYStressGP[iTimeGP][iBndGP] - impAndEta[ltsFace].eta_s * slipRateStrike[ltsFace][iBndGP];
        faultStresses.XZTractionResultGP[iTimeGP][iBndGP] = faultStresses.XZStressGP[iTimeGP][iBndGP] - impAndEta[ltsFace].eta_s * slipRateDip[ltsFace][iBndGP];
        tractionXY[ltsFace][iBndGP] = faultStresses.XYTractionResultGP[iTimeGP][iBndGP];
        tractionXZ[ltsFace][iBndGP] = faultStresses.XYTractionResultGP[iTimeGP][iBndGP];

        //-------------------------------------
        //update Directional Slip
        slipStrike[ltsFace][iBndGP] += slipRateStrike[ltsFace][iBndGP] * deltaT[iTimeGP];
        slipDip[ltsFace][iBndGP] += slipRateDip[ltsFace][iBndGP] * deltaT[iTimeGP];
      }
    }


/*
 * evaluate friction law: updated mu -> friction law
 * for example see Carsten Uphoff's thesis: Eq. 2.45
 */
    virtual void frictionFunctionHook(std::array<real, numOfPointsPadded> &stateVariablePsi, unsigned int ltsFace) {
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
        mu[ltsFace][iBndGP] = mu_S[ltsFace][iBndGP] - (mu_S[ltsFace][iBndGP] - mu_D[ltsFace][iBndGP]) * stateVariablePsi[iBndGP];
      }
    }


    /*
     * instantaneous healing option Reset Mu and Slip
     */
    virtual void instantaneousHealing(unsigned int ltsFace){
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
        if (SlipRateMagnitude[ltsFace][iBndGP] < u_0) {
          mu[ltsFace][iBndGP] = mu_S[ltsFace][iBndGP];
          slip[ltsFace][iBndGP] = 0.0;
        }
      }
    }
/*
 * output time when shear stress is equal to the dynamic stress after rupture arrived
 * currently only for linear slip weakening
 */
  virtual void saveDynamicStressOutput(
      unsigned int ltsFace
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      if (rupture_time[ltsFace][iBndGP] > 0.0 &&
          rupture_time[ltsFace][iBndGP] <= m_fullUpdateTime &&
          DS[iBndGP] &&
          std::fabs(slip[ltsFace][iBndGP]) >= d_c[ltsFace][iBndGP]) {
        dynStress_time[ltsFace][iBndGP] = m_fullUpdateTime;
        DS[ltsFace][iBndGP] = false;
      }
    }
  }


};//End of Class LinearSlipWeakeningLaw

class seissol::dr::friction_law::LinearSlipWeakeningLawFL2 : public seissol::dr::friction_law::LinearSlipWeakeningLaw<seissol::dr::friction_law::LinearSlipWeakeningLawFL2> {
public:
  /*
   * compute fault strength
   */
  virtual void calcStrengthHook(std::array<real, numOfPointsPadded> &Strength,
                                FaultStresses &faultStresses,
                                unsigned int iTimeGP,
                                unsigned int ltsFace);

  /*
   * compute state variable
   */
  virtual void calcStateVariableHook(std::array<real, numOfPointsPadded> &stateVariablePsi,
                                     std::array<real, numOfPointsPadded> &outputSlip,
                                     dynamicRupture::kernel::resampleParameter &resampleKrnl,
                                     unsigned int iTimeGP,
                                     unsigned int ltsFace);
};//End of Class LinearSlipWeakeningLawFL2

class seissol::dr::friction_law::LinearSlipWeakeningLawFL16 : public seissol::dr::friction_law::LinearSlipWeakeningLawFL2 {
protected:
  real                    (*forced_rupture_time)[numOfPointsPadded];
  real                    *tn;

  /*
* copies all parameters from the DynamicRupture LTS to the local attributes
*/
  void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) override;

  virtual void setTimeHook(unsigned int ltsFace) override;

  virtual void calcStateVariableHook(std::array<real, numOfPointsPadded> &stateVariablePsi,
                                     std::array<real, numOfPointsPadded> &outputSlip,
                                     dynamicRupture::kernel::resampleParameter &resampleKrnl,
                                     unsigned int iTimeGP,
                                     unsigned int ltsFace) override;
}; //End of Class LinearSlipWeakeningLawFL16

/*
 * Law for Bimaterial faults, implements strength regularization (according to prakash clifton)
 * currently regularized strength is not used (bug)
 * State variable (slip) is not resampled in this friction law!
 */
class seissol::dr::friction_law::LinearSlipWeakeningLawBimaterialFL6 : public seissol::dr::friction_law::LinearSlipWeakeningLaw<seissol::dr::friction_law::LinearSlipWeakeningLawBimaterialFL6> {
public:
  virtual void calcStrengthHook(std::array<real, numOfPointsPadded> &Strength,
      FaultStresses &faultStresses,
      unsigned int iTimeGP,
      unsigned int ltsFace);

  virtual void calcStateVariableHook(std::array<real, numOfPointsPadded> &stateVariablePsi,
      std::array<real, numOfPointsPadded> &outputSlip,
      dynamicRupture::kernel::resampleParameter &resampleKrnl,
      unsigned int iTimeGP,
      unsigned int ltsFace);

protected:
  //Attributes
  real  (*strengthData)[numOfPointsPadded];

  /*
 * copies all parameters from the DynamicRupture LTS to the local attributes
 */
  void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) override;

  /*
   * calculates strength
   */
  void prak_clif_mod(real &strength, real &sigma, real &LocSlipRate, real &mu, real &dt);
};


#endif //SEISSOL_LINEARSLIPWEAKENING_H
