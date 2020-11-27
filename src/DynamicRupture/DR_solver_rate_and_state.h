//
// Created by adrian on 03.09.20.
//

#ifndef SEISSOL_DR_SOLVER_RATE_AND_STATE_H
#define SEISSOL_DR_SOLVER_RATE_AND_STATE_H

#include "DR_solver_base.h"

#include "DR_math.h"
#include <yaml-cpp/yaml.h>



namespace seissol {
  namespace dr {
    namespace fr_law {
      template <class Derived>
      class RateAndStateSolver;
      class RateAndStateNucFL103;  //rate and state slip law with nucleation
      class RateAndStateThermalFL103;
    }
  }
}

template <class Derived>
class seissol::dr::fr_law::RateAndStateSolver : public seissol::dr::fr_law::BaseFrictionSolver {

protected:
  //!PARAMETERS of THE optimisation loops
  //!absolute tolerance on the function to be optimzed
  //! This value is quite arbitrary (a bit bigger as the expected numerical error) and may not be the most adapted
  //! Number of iteration in the loops
  const unsigned int nSRupdates = 60;
  const unsigned int nSVupdates = 2;



public:
  virtual void evaluate(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real fullUpdateTime,
                        real timeWeights[CONVERGENCE_ORDER]) override {
    //first copy all Variables from the Base Lts dynRup tree
    static_cast<Derived*>(this)->copyLtsTreeToLocalRS(layerData, dynRup, fullUpdateTime);

    static_cast<Derived*>(this)->preCalcTime();


#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {

      //initialize local variables inside parallel face loop
      bool has_converged = false;
      FaultStresses faultStresses = {};
      real deltaStateVar[numOfPointsPadded] = {0};
      std::array<real, numOfPointsPadded> tmpSlip{0};   //required for averageSlip calculation
      std::array<real, numOfPointsPadded> normalStress{0};
      std::array<real, numOfPointsPadded> TotalShearStressYZ{0};
      std::array<real, numOfPointsPadded> stateVarZero{0};
      std::array<real, numOfPointsPadded> SR_tmp{0};
      std::array<real, numOfPointsPadded> LocSV{0};
      std::array<real, numOfPointsPadded> SRtest{0};

      //for thermalPressure
      std::array<real, numOfPointsPadded> P_f{0};

      precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      static_cast<Derived*>(this)->setInitialValues(LocSV, ltsFace);

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {

        static_cast<Derived*>(this)->calcInitialSlipRate(TotalShearStressYZ,faultStresses, stateVarZero, LocSV, SR_tmp, iTimeGP, ltsFace);
        static_cast<Derived*>(this)->hookSetInitialP_f(P_f, ltsFace);

        for (unsigned int j = 0; j < nSVupdates; j++) {
          static_cast<Derived*>(this)->hookCalcP_f(P_f, faultStresses, false, iTimeGP, ltsFace);
          static_cast<Derived*>(this)->updateStateVariableIterative(has_converged, stateVarZero, SR_tmp, LocSV, P_f,
              normalStress, TotalShearStressYZ, SRtest, faultStresses, iTimeGP, ltsFace);
        } //End nSVupdates-loop   j=1,nSVupdates   !This loop corrects SV values
        if (!has_converged) {
          static_cast<Derived*>(this)->executeIfNotConverged(LocSV, ltsFace);
        }
        static_cast<Derived*>(this)->hookCalcP_f(P_f, faultStresses, true, iTimeGP, ltsFace);
        static_cast<Derived*>(this)->calcSlipRateAndTraction(stateVarZero, SR_tmp, LocSV, normalStress,
            TotalShearStressYZ, tmpSlip, deltaStateVar, faultStresses, iTimeGP, ltsFace);

      } // End of iTimeGP-loop

      static_cast<Derived*>(this)->resampleStateVar(deltaStateVar, ltsFace);

      //---------------------------------------------

      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      saveRuptureFrontOutput(ltsFace);

      savePeakSlipRateOutput(ltsFace);

      //output time when shear stress is equal to the dynamic stress after rupture arrived
      //currently only for linear slip weakening
      static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace);

      //---compute and store slip to determine the magnitude of an earthquake ---
      //    to this end, here the slip is computed and averaged per element
      //    in calc_seissol.f90 this value will be multiplied by the element surface
      //    and an output happened once at the end of the simulation
      saveAverageSlipOutput(tmpSlip, ltsFace);

      postcomputeImposedStateFromNewStress(
          QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace],
          faultStresses, timeWeights, ltsFace);

    }//end face loop
  }//end evaluate function
};


class seissol::dr::fr_law::RateAndStateNucFL103 : public seissol::dr::fr_law::RateAndStateSolver<seissol::dr::fr_law::RateAndStateNucFL103> {
protected:
  //Attributes
  real  (*nucleationStressInFaultCS)[numOfPointsPadded][6];
  real dt = 0;
  real Gnuc = 0;

  real  (*RS_a_array)[numOfPointsPadded];
  real  (*RS_srW_array)[numOfPointsPadded];
  real  (*RS_sl0_array)[numOfPointsPadded];

  bool  (*DS)[numOfPointsPadded];
  real  (*stateVar)[numOfPointsPadded];
  real  (*dynStress_time)[numOfPointsPadded];

  //!TU 7.07.16: if the SR is too close to zero, we will have problems (NaN)
  //!as a consequence, the SR is affected the AlmostZero value when too small
  const double AlmostZero = 1e-45;

  //!PARAMETERS of THE optimisation loops
  //!absolute tolerance on the function to be optimzed
  //! This value is quite arbitrary (a bit bigger as the expected numerical error) and may not be the most adapted
  //! Number of iteration in the loops
  const unsigned int nSRupdates = 60;
  const unsigned int nSVupdates = 2;

  const double aTolF = 1e-8;

public:
  /*
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  virtual void copyLtsTreeToLocalRS(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) {
    //first copy all Variables from the Base Lts dynRup tree
    BaseFrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    //maybe change later to const_cast?
    seissol::initializers::DR_FL_103 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103 *>(dynRup);
    nucleationStressInFaultCS =  layerData.var(ConcreteLts->nucleationStressInFaultCS); ;

    RS_sl0_array    = layerData.var(ConcreteLts->RS_sl0_array);
    RS_a_array      = layerData.var(ConcreteLts->RS_a_array);
    RS_srW_array    = layerData.var(ConcreteLts->RS_srW_array);
    DS              = layerData.var(ConcreteLts->DS);
    averaged_Slip   = layerData.var(ConcreteLts->averaged_Slip);
    stateVar        = layerData.var(ConcreteLts->stateVar);
    dynStress_time  = layerData.var(ConcreteLts->dynStress_time);
  }

  void preCalcTime(){
    dt = 0;
    for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {
      dt += deltaT[iTimeGP];
    }
    if (m_fullUpdateTime <= m_Params->t_0) {
      Gnuc = Calc_SmoothStepIncrement(m_fullUpdateTime, dt);
    }
  }

  void setInitialValues(std::array<real, numOfPointsPadded> &LocSV, unsigned int ltsFace){
    if (m_fullUpdateTime <= m_Params->t_0) {
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
        for (int i = 0; i < 6; i++) {
          initialStressInFaultCS[ltsFace][iBndGP][i] += nucleationStressInFaultCS[ltsFace][iBndGP][i] * Gnuc;
        }
      }
    } //end If-Tnuc


    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      LocSV[iBndGP] = stateVar[ltsFace][iBndGP];     //DISC%DynRup%StateVar(iBndGP,iFace)      //local varriable required
    }
  }

  void calcInitialSlipRate(
      std::array<real, numOfPointsPadded> &TotalShearStressYZ,
      FaultStresses &faultStresses,
      std::array<real, numOfPointsPadded> &stateVarZero,
      std::array<real, numOfPointsPadded> &LocSV,
      std::array<real, numOfPointsPadded> &SR_tmp,
      unsigned int iTimeGP,
      unsigned int ltsFace){

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      // friction develops as                    mu = a * arcsinh[ V/(2*V0) * exp(SV/a) ]
      // state variable SV develops as          dSV / dt = -(V - L) * (SV - SV_ss)
      //                                        SV_ss = a * ln[ 2*V0/V * sinh(mu_ss/a) ]
      //                                        mu_ss = mu_w + [mu_lv - mu_w] / [ 1 + (V/Vw)^8 ] ^ (1/8) ]
      //                                        mu_lv = mu_0 - (b-a) ln (V/V0)

      TotalShearStressYZ[iBndGP] = std::sqrt(
          seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP], 2) +
          seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP], 2));

      // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
      // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )
      stateVarZero[iBndGP] = LocSV[iBndGP];    // Careful, the SV must always be corrected using SV0 and not LocSV!

      // The following process is adapted from that described by Kaneko et al. (2008)
      SlipRateMagnitude[ltsFace][iBndGP] = std::sqrt(seissol::dr::aux::power(slipRateStrike[ltsFace][iBndGP], 2) + seissol::dr::aux::power(slipRateDip[ltsFace][iBndGP], 2) );
      SlipRateMagnitude[ltsFace][iBndGP] = std::max(AlmostZero, SlipRateMagnitude[ltsFace][iBndGP]);
      SR_tmp[iBndGP] = SlipRateMagnitude[ltsFace][iBndGP];

    }// End of iBndGP-loop
  }


  void updateStateVariableIterative(bool &has_converged,
                                    std::array<real, numOfPointsPadded> &stateVarZero,
                                    std::array<real, numOfPointsPadded> &SR_tmp,
                                    std::array<real, numOfPointsPadded> &LocSV,
                                    std::array<real, numOfPointsPadded> &P_f,
                                    std::array<real, numOfPointsPadded> &normalStress,
                                    std::array<real, numOfPointsPadded> &TotalShearStressYZ,
                                    std::array<real, numOfPointsPadded> &SRtest,
                                    FaultStresses &faultStresses,
                                    unsigned int iTimeGP,
                                    unsigned int ltsFace){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      //fault strength using LocMu and P_f from previous timestep/iteration
      //1.update SV using Vold from the previous time step
      updateStateVariable(iBndGP, ltsFace, stateVarZero[iBndGP], deltaT[iTimeGP], SR_tmp[iBndGP], LocSV[iBndGP]);
      normalStress[iBndGP] = faultStresses.NormalStressGP[iTimeGP][iBndGP] + initialStressInFaultCS[ltsFace][iBndGP][0] - P_f[iBndGP];
    }// End of iBndGP-loop

    //2. solve for Vnew , applying the Newton-Raphson algorithm
    //effective normal stress including initial stresses and pore fluid pressure
    has_converged = IterativelyInvertSR(ltsFace, nSRupdates, LocSV, normalStress, TotalShearStressYZ, SRtest);

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      // 3. update theta, now using V=(Vnew+Vold)/2
      // For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
      SR_tmp[iBndGP] = 0.5 * (SlipRateMagnitude[ltsFace][iBndGP] + fabs(SRtest[iBndGP]));

      // 4. solve again for Vnew
      SlipRateMagnitude[ltsFace][iBndGP] = fabs(SRtest[iBndGP]);

      //!update LocMu
      updateMu(ltsFace, iBndGP, LocSV[iBndGP]);
    }// End of iBndGP-loop
  }

  void executeIfNotConverged(std::array<real, numOfPointsPadded> &LocSV, unsigned ltsFace){
    real tmp = 0.5 / m_Params->rs_sr0 * exp(LocSV[0] / RS_a_array[ltsFace][0]) * SlipRateMagnitude[ltsFace][0];
    //!logError(*) 'nonConvergence RS Newton', time
    std::cout << "nonConvergence RS Newton, time: " << m_fullUpdateTime << std::endl;
    assert(!std::isnan(tmp) && "nonConvergence RS Newton");
  }

  void calcSlipRateAndTraction(
      std::array<real, numOfPointsPadded> &stateVarZero,
      std::array<real, numOfPointsPadded> &SR_tmp,
      std::array<real, numOfPointsPadded> &LocSV,
      std::array<real, numOfPointsPadded> &normalStress,
      std::array<real, numOfPointsPadded> &TotalShearStressYZ,
      std::array<real, numOfPointsPadded> &tmpSlip,
      real deltaStateVar[numOfPointsPadded],
      FaultStresses &faultStresses,
      unsigned int iTimeGP,
      unsigned int ltsFace){
    std::array<real, numOfPointsPadded> LocSlipRateMagnitude{0};

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      //! SV from mean slip rate in tmp
      updateStateVariable(iBndGP, ltsFace, stateVarZero[iBndGP], deltaT[iTimeGP], SR_tmp[iBndGP], LocSV[iBndGP]);

      //!update LocMu for next strength determination, only needed for last update
      updateMu(ltsFace, iBndGP, LocSV[iBndGP]);

      //! update stress change
      tracXY[ltsFace][iBndGP] = -((initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP]) / TotalShearStressYZ[iBndGP]) * mu[ltsFace][iBndGP] * normalStress[iBndGP];
      tracXZ[ltsFace][iBndGP] = -((initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP]) / TotalShearStressYZ[iBndGP]) * mu[ltsFace][iBndGP] * normalStress[iBndGP];
      tracXY[ltsFace][iBndGP] -= initialStressInFaultCS[ltsFace][iBndGP][3];
      tracXZ[ltsFace][iBndGP] -= initialStressInFaultCS[ltsFace][iBndGP][5];

      //Compute slip
      //! ABS of locSlipRate removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening
      slip[ltsFace][iBndGP] += SlipRateMagnitude[ltsFace][iBndGP] * deltaT[iTimeGP];

      //!Update slip rate (notice that locSlipRate(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
      slipRateStrike[ltsFace][iBndGP] = -impAndEta[ltsFace].inv_eta_s * (tracXY[ltsFace][iBndGP] - faultStresses.XYStressGP[iTimeGP][iBndGP]);
      slipRateDip[ltsFace][iBndGP] = -impAndEta[ltsFace].inv_eta_s * (tracXZ[ltsFace][iBndGP] - faultStresses.XZStressGP[iTimeGP][iBndGP]);

      //!TU 07.07.16: correct slipRateStrike and slipRateDip to avoid numerical errors
      LocSlipRateMagnitude[iBndGP] = sqrt(seissol::dr::aux::power(slipRateStrike[ltsFace][iBndGP], 2) + seissol::dr::aux::power(slipRateDip[ltsFace][iBndGP], 2));
      if (LocSlipRateMagnitude[iBndGP] != 0) {
        slipRateStrike[ltsFace][iBndGP] *= SlipRateMagnitude[ltsFace][iBndGP] / LocSlipRateMagnitude[iBndGP];
        slipRateDip[ltsFace][iBndGP] *= SlipRateMagnitude[ltsFace][iBndGP] / LocSlipRateMagnitude[iBndGP];
      }
      tmpSlip[iBndGP] += LocSlipRateMagnitude[iBndGP] * deltaT[iTimeGP];

      slip1[ltsFace][iBndGP] += slipRateStrike[ltsFace][iBndGP] * deltaT[iTimeGP];
      slip2[ltsFace][iBndGP] += slipRateDip[ltsFace][iBndGP] * deltaT[iTimeGP];

      //!Save traction for flux computation
      faultStresses.XYTractionResultGP[iTimeGP][iBndGP] = tracXY[ltsFace][iBndGP];
      faultStresses.XZTractionResultGP[iTimeGP][iBndGP] = tracXZ[ltsFace][iBndGP];

      //Could be outside TimeLoop, since only last time result is used later
      deltaStateVar[iBndGP] = LocSV[iBndGP] - stateVar[ltsFace][iBndGP];
    } // End of BndGP-loop

  }

  void resampleStateVar(
      real deltaStateVar[numOfPointsPadded],
      unsigned int ltsFace){

    dynamicRupture::kernel::resampleParameter resampleKrnl = {};
    resampleKrnl.resampleM = init::resample::Values;
    real resampledDeltaStateVar[numOfPointsPadded] = {0};
    resampleKrnl.resamplePar = deltaStateVar;
    resampleKrnl.resampledPar = resampledDeltaStateVar;  //output from execute
    resampleKrnl.execute();

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      //write back State Variable to lts tree
      stateVar[ltsFace][iBndGP] = stateVar[ltsFace][iBndGP] + resampledDeltaStateVar[iBndGP];
    }
  }


  //output time when shear stress is equal to the dynamic stress after rupture arrived
  //currently only for linear slip weakening
  void saveDynamicStressOutput(
      unsigned int face
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      if (rupture_time[face][iBndGP] > 0.0 &&
          rupture_time[face][iBndGP] <= m_fullUpdateTime &&
          DS[iBndGP] &&
          mu[face][iBndGP] <= ( m_Params->mu_w+0.05*(m_Params->rs_f0-m_Params->mu_w) ) ) {
        dynStress_time[face][iBndGP] = m_fullUpdateTime;
        DS[face][iBndGP] = false;
      }
    }
  }

  virtual void hookSetInitialP_f(std::array<real, numOfPointsPadded> &P_f, unsigned int ltsFace){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      P_f[iBndGP] = 0.0;
    }
  }

  virtual void hookCalcP_f(std::array<real, numOfPointsPadded> &P_f, FaultStresses &faultStresses, bool saveTmpInTP, unsigned int iTimeGP, unsigned int ltsFace){
  }

protected:
  void updateStateVariable(int iBndGP, unsigned int face, real SV0, real time_inc, real &SR_tmp, real &LocSV){
    double flv, fss, SVss;
    double RS_fw = m_Params->mu_w;
    double RS_srW = RS_srW_array[face][iBndGP];
    double RS_a = RS_a_array[face][iBndGP];
    double RS_sl0 = RS_sl0_array[face][iBndGP];
    double exp1;

    // low-velocity steady state friction coefficient
    flv = m_Params->rs_f0 - (m_Params->rs_b-RS_a)* log(SR_tmp/m_Params->rs_sr0);
    // steady state friction coefficient
    fss = RS_fw + (flv - RS_fw)/pow(1.0+seissol::dr::aux::power(SR_tmp/RS_srW,8.0) ,1.0/8.0);
    // steady-state state variable
    // For compiling reasons we write SINH(X)=(EXP(X)-EXP(-X))/2
    SVss = RS_a * log(2.0*m_Params->rs_sr0/SR_tmp * (exp(fss/RS_a)-exp(-fss/RS_a))/2.0);

    // exact integration of dSV/dt DGL, assuming constant V over integration step

    exp1 = exp(-SR_tmp*(time_inc/RS_sl0) );
    LocSV = SVss*(1.0-exp1)+exp1*SV0;

    assert( !(std::isnan(LocSV) && iBndGP < numberOfPoints ) && "NaN detected");
  }

  /*
   * If the function did not converge it returns false
   * Newton method to solve non-linear equation of slip rate
   * solution returned in SRtest
   */
  bool IterativelyInvertSR (unsigned int ltsFace, int nSRupdates,
                            std::array<real, numOfPointsPadded> &LocSV, std::array<real, numOfPointsPadded> &n_stress,
                            std::array<real, numOfPointsPadded> &sh_stress, std::array<real, numOfPointsPadded> &SRtest ){

    double tmp[numOfPointsPadded], tmp2[numOfPointsPadded], tmp3[numOfPointsPadded], mu_f[numOfPointsPadded], dmu_f[numOfPointsPadded], NR[numOfPointsPadded], dNR[numOfPointsPadded];
    //double AlmostZero = 1e-45;
    bool has_converged = false;

    //!solve for Vnew = SR , applying the Newton-Raphson algorithm
    //!SR fulfills g(SR)=f(SR)
    //!-> find root of NR=f-g using a Newton-Raphson algorithm with dNR = d(NR)/d(SR)
    //!SR_{i+1}=SR_i-( NR_i / dNR_i )
    //!
    //!        equalize:
    //!         g = SR*MU/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
    //!         f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
    //!  where mu = friction coefficient, dependening on the RSF law used

    for(int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++){
      //! first guess = SR value of the previous step
      SRtest[iBndGP] = SlipRateMagnitude[ltsFace][iBndGP];
      tmp[iBndGP]   =  0.5 / m_Params->rs_sr0 *exp(LocSV[iBndGP]/RS_a_array[ltsFace][iBndGP]);
    }

    for(int i = 0; i < nSRupdates; i++){
      for(int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++){

        //!f = ( tmp2 * ABS(LocP+P_0)- ABS(S_0))*(S_0)/ABS(S_0)
        //!g = SRtest * 1.0/(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) + ABS(ShTest)
        //!for compiling reasons ASINH(X)=LOG(X+SQRT(X^2+1))

        //!calculate friction coefficient
        tmp2[iBndGP]  = tmp[iBndGP]*SRtest[iBndGP];
        mu_f[iBndGP]  = RS_a_array[ltsFace][iBndGP] * log(tmp2[iBndGP] + sqrt(seissol::dr::aux::power(tmp2[iBndGP], 2) + 1.0));
        dmu_f[iBndGP] = RS_a_array[ltsFace][iBndGP] / sqrt(1.0 + seissol::dr::aux::power(tmp2[iBndGP], 2)) * tmp[iBndGP];
        NR[iBndGP]    = -impAndEta[ltsFace].inv_eta_s * (fabs(n_stress[iBndGP]) * mu_f[iBndGP] - sh_stress[iBndGP]) - SRtest[iBndGP];
      }

      has_converged = true;

      //max element of NR must be smaller then aTolF
      for(int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++){
        if (fabs(NR[iBndGP]) >= aTolF ){
          has_converged = false;
          break;
        }
      }
      if(has_converged){
        return has_converged;
      }
      for(int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++){

        //!derivative of NR
        dNR[iBndGP]   = -impAndEta[ltsFace].inv_eta_s * (fabs(n_stress[iBndGP]) * dmu_f[iBndGP]) - 1.0;
        //!ratio
        tmp3[iBndGP] = NR[iBndGP]/dNR[iBndGP];

        //!update SRtest
        SRtest[iBndGP] = std::max(AlmostZero,SRtest[iBndGP]-tmp3[iBndGP]);
      }
    }
  }

  /*
   * Alternative to IterativelyInvertSR: Instead of newton, brents method is used:
   * slower convergence but higher stability
   * if the function did not converge it returns false
   * From Upoffc git Zero.ccp (https://github.com/TEAR-ERC/tandem/blob/master/src/util/Zero.cpp)
 */
  bool IterativelyInvertSR_Brent(unsigned int ltsFace, int nSRupdates,
                                 std::array<real, numOfPointsPadded> &LocSV, std::array<real, numOfPointsPadded> &n_stress,
                                 std::array<real, numOfPointsPadded> &sh_stress,  std::array<real, numOfPointsPadded> &SRtest ){
    std::function<double(double, int)> F;
    double tol = 1e-30;

    double *RS_a = RS_a_array[ltsFace];
    double RS_sr0_ = m_Params->rs_sr0;
    double invEta = impAndEta[ltsFace].inv_eta_s;

    F = [invEta, &sh_stress, n_stress, RS_a, LocSV, RS_sr0_](double SR, int iBndGP){
      double tmp  = 0.5 / RS_sr0_ *exp(LocSV[iBndGP]/RS_a[iBndGP]) * SR;
      double mu_f  = RS_a[iBndGP] * log(tmp+sqrt(seissol::dr::aux::power(tmp,2)+1.0));
      return -invEta * (fabs(n_stress[iBndGP])*mu_f-sh_stress[iBndGP])-SR;
    };

    for(int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++){
      //TODO: change boundaries?
      double a = SlipRateMagnitude[ltsFace][iBndGP] - impAndEta[ltsFace].inv_eta_s * sh_stress[iBndGP];
      double b  = SlipRateMagnitude[ltsFace][iBndGP] + impAndEta[ltsFace].inv_eta_s * sh_stress[iBndGP];

      double eps = std::numeric_limits<double>::epsilon();
      double Fa = F(a, iBndGP);
      //if(std::isinf(Fa)){
      //  Fa = std::numeric_limits<double>::max();
      //}
      double Fb = F(b, iBndGP);
      assert(std::copysign(Fa, Fb) != Fa); // Fa and Fb have different signs
      double c = a;
      double Fc = Fa;
      double d = b - a;
      double e = d;
      while (Fb != 0.0) {
        if (std::copysign(Fb, Fc) == Fb) {
          c = a;
          Fc = Fa;
          d = b - a;
          e = d;
        }
        if (std::fabs(Fc) < std::fabs(Fb)) {
          a = b;
          b = c;
          c = a;
          Fa = Fb;
          Fb = Fc;
          Fc = Fa;
        }
        // Convergence test
        double xm = 0.5 * (c - b);
        double tol1 = 2.0 * eps * std::fabs(b) + 0.5 * tol;
        if (std::fabs(xm) <= tol1 || Fb == 0.0) {
          break;
        }
        if (std::fabs(e) < tol1 || std::fabs(Fa) <= std::fabs(Fb)) {
          // bisection
          d = xm;
          e = d;
        } else {
          double s = Fb / Fa;
          double p, q;
          if (a != c) {
            // linear interpolation
            q = Fa / Fc;
            double r = Fb / Fc;
            p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
            q = (q - 1.0) * (r - 1.0) * (s - 1.0);
          } else {
            // inverse quadratic interpolation
            p = 2.0 * xm * s;
            q = 1.0 - s;
          }
          if (p > 0) {
            q = -q;
          } else {
            p = -p;
          }
          if (2.0 * p < 3.0 * xm * q - std::fabs(tol1 * q) && p < std::fabs(0.5 * e * q)) {
            e = d;
            d = p / q;
          } else {
            // bisection
            d = xm;
            e = d;
          }
        }
        a = b;
        Fa = Fb;
        if (std::fabs(d) > tol1) {
          b += d;
        } else {
          b += std::copysign(tol1, xm);
        }
        Fb = F(b, iBndGP);
      }
      SRtest[iBndGP] = b;
    }
    return true;
  }
  /*
   * mu = a * arcsinh[ V/(2*V0) * exp(SV/a) ]
   */
  void updateMu(unsigned int ltsFace, unsigned int iBndGP, real LocSV){
    //! X in Asinh(x) for mu calculation
    real tmp = 0.5 / m_Params->rs_sr0 * exp(LocSV / RS_a_array[ltsFace][iBndGP]) * SlipRateMagnitude[ltsFace][iBndGP];
    //! mu from locSlipRate with SINH(X)=LOG(X+SQRT(X^2+1))
    mu[ltsFace][iBndGP] = RS_a_array[ltsFace][iBndGP] * log(tmp + sqrt(seissol::dr::aux::power(tmp, 2) + 1.0));
  }
}; //end class FL_103


class seissol::dr::fr_law::RateAndStateThermalFL103 : public seissol::dr::fr_law::RateAndStateNucFL103 {
protected:

  real (*temperature)[numOfPointsPadded];
  real (*pressure)[numOfPointsPadded];
  real (*TP_Theta)[numOfPointsPadded][TP_grid_nz];
  real (*TP_sigma)[numOfPointsPadded][TP_grid_nz];
  real (*TP_half_width_shear_zone)[numOfPointsPadded];
  real (*alpha_hy)[numOfPointsPadded];

  real TP_grid[TP_grid_nz];
  real TP_DFinv[TP_grid_nz];

  real Sh[numOfPointsPadded];
  real Theta_tmp[TP_grid_nz];
  real Sigma_tmp[TP_grid_nz];

public:
  void initializeTP(seissol::Interoperability &e_interoperability){
    e_interoperability.getDynRupTP(TP_grid, TP_DFinv);
  }

  /*
 * copies all parameters from the DynamicRupture LTS to the local attributes
 */
  virtual void copyLtsTreeToLocalRS(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) override {
    //first copy all Variables from the Base Lts dynRup tree
    RateAndStateNucFL103::copyLtsTreeToLocalRS(layerData, dynRup, fullUpdateTime);

    //maybe change later to const_cast?
    seissol::initializers::DR_FL_103_Thermal *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103_Thermal *>(dynRup);
    temperature               = layerData.var(ConcreteLts->temperature);
    pressure                  = layerData.var(ConcreteLts->pressure);
    TP_Theta                  = layerData.var(ConcreteLts->TP_theta);
    TP_sigma                  = layerData.var(ConcreteLts->TP_sigma);
    TP_half_width_shear_zone  = layerData.var(ConcreteLts->TP_half_width_shear_zone);
    alpha_hy                  = layerData.var(ConcreteLts->alpha_hy);
  }

protected:
  void hookSetInitialP_f(std::array<real, numOfPointsPadded> &P_f, unsigned int ltsFace) override{
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      P_f[iBndGP] = pressure[ltsFace][iBndGP];
    }
  }

  void hookCalcP_f(std::array<real, numOfPointsPadded> &P_f,  FaultStresses &faultStresses, bool saveTmpInTP, unsigned int iTimeGP, unsigned int ltsFace) override {
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      Sh[iBndGP] = -mu[ltsFace][iBndGP] * (faultStresses.NormalStressGP[iTimeGP][iBndGP] + initialStressInFaultCS[ltsFace][iBndGP][0] - P_f[iBndGP]);

      for (unsigned int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++) {
        //!recover original values as it gets overwritten in the ThermalPressure routine
        Theta_tmp[iTP_grid_nz] = TP_Theta[ltsFace][iBndGP][iTP_grid_nz];
        Sigma_tmp[iTP_grid_nz] = TP_sigma[ltsFace][iBndGP][iTP_grid_nz];
      }
      //!use Theta/Sigma from last call in this update, dt/2 and new SR from NS
      Calc_ThermalPressure(iBndGP, iTimeGP, ltsFace);

      P_f[iBndGP] = pressure[ltsFace][iBndGP];
      if(saveTmpInTP){
        for (unsigned int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++) {
          TP_Theta[ltsFace][iBndGP][iTP_grid_nz] = Theta_tmp[iTP_grid_nz];
          TP_sigma[ltsFace][iBndGP][iTP_grid_nz] = Sigma_tmp[iTP_grid_nz];
        }
      }
    }
  }

  void Calc_ThermalPressure(unsigned int iBndGP, unsigned int iTimeGP, unsigned int ltsFace){
    real tauV, Lambda_prime, T, p;
    real tmp[TP_grid_nz];
    real omega_theta[TP_grid_nz];
    real omega_sigma[TP_grid_nz];
    real theta_current[TP_grid_nz];
    real sigma_current[TP_grid_nz];

    T = 0.0;
    p = 0.0;

    tauV = Sh[iBndGP] * SlipRateMagnitude[ltsFace][iBndGP]; //!fault strenght*slip rate
    Lambda_prime = m_Params->TP_lambda * m_Params->alpha_th / (alpha_hy[ltsFace][iBndGP] - m_Params->alpha_th);

    for (unsigned int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++) {
      //!Gaussian shear zone in spectral domain, normalized by w
      tmp[iTP_grid_nz] = seissol::dr::aux::power(TP_grid[iTP_grid_nz] / TP_half_width_shear_zone[ltsFace][iBndGP], 2);
      //!1. Calculate diffusion of the field at previous timestep

      //!temperature
      theta_current[iTP_grid_nz] = Theta_tmp[iTP_grid_nz] * exp(-m_Params->alpha_th * deltaT[iTimeGP] * tmp[iTP_grid_nz]);
      //!pore pressure + lambda'*temp
      sigma_current[iTP_grid_nz] = Sigma_tmp[iTP_grid_nz] * exp(-alpha_hy[ltsFace][iBndGP] * deltaT[iTimeGP] * tmp[iTP_grid_nz]);

      //!2. Add current contribution and get new temperature
      omega_theta[iTP_grid_nz] = heat_source(tmp[iTP_grid_nz], m_Params->alpha_th, iTP_grid_nz, iTimeGP);
      Theta_tmp[iTP_grid_nz] = theta_current[iTP_grid_nz] + (tauV / m_Params->rho_c) * omega_theta[iTP_grid_nz];
      omega_sigma[iTP_grid_nz] = heat_source(tmp[iTP_grid_nz], alpha_hy[ltsFace][iBndGP], iTP_grid_nz, iTimeGP);
      Sigma_tmp[iTP_grid_nz] = sigma_current[iTP_grid_nz] + ((m_Params->TP_lambda + Lambda_prime) * tauV) / (m_Params->rho_c) * omega_sigma[iTP_grid_nz];

      //!3. Recover temperature and pressure using inverse Fourier
      //! transformation with the calculated fourier coefficients

      //!new contribution
      T += (TP_DFinv[iTP_grid_nz] / TP_half_width_shear_zone[ltsFace][iBndGP]) * Theta_tmp[iTP_grid_nz];
      p += (TP_DFinv[iTP_grid_nz] / TP_half_width_shear_zone[ltsFace][iBndGP]) * Sigma_tmp[iTP_grid_nz];
    }
    //Update pore pressure change (sigma = pore pressure + lambda'*temp)
    //In the BIEM code (Lapusta) they use T without initial value
    p = p - Lambda_prime*T;

    //Temp and pore pressure change at single GP on the fault + initial values
    temperature[ltsFace][iBndGP] = T + m_Params->IniTemp;
    pressure[ltsFace][iBndGP] = -p + m_Params->IniPressure;
  }

  real heat_source(real tmp, real alpha, unsigned int iTP_grid_nz, unsigned int iTimeGP){
    //!original function in spatial domain
    //!omega = 1/(w*sqrt(2*pi))*exp(-0.5*(z/TP_half_width_shear_zone).^2);
    //!function in the wavenumber domain *including additional factors in front of the heat source function*
    //!omega = 1/(*alpha*Dwn**2**(sqrt(2.0*pi))*exp(-0.5*(Dwn*TP_half_width_shear_zone)**2)*(1-exp(-alpha**dt**tmp))
    //!inserting Dwn/TP_half_width_shear_zone (scaled) for Dwn cancels out TP_half_width_shear_zone
    return 1.0/(alpha*tmp*(sqrt(2.0*M_PI))) * exp(-0.5*seissol::dr::aux::power(TP_grid[iTP_grid_nz], 2)) * (1.0 - exp(-alpha*deltaT[iTimeGP]*tmp));
  }
};


#endif //SEISSOL_DR_SOLVER_RATE_AND_STATE_H
