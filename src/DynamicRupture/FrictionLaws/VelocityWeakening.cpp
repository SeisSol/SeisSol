//
// Created by sebastian on 13/10/2021.
//

#include "VelocityWeakening.h"

void seissol::dr::friction_law::VelocityWeakening::copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                                                                      seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) {
  //first copy all Variables from the Base Lts dynRup tree
  BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  //TODO: change later to const_cast
  //seissol::initializers::DR_lts_template *ConcreteLts = dynamic_cast<seissol::initializers::DR_lts_template *>(dynRup);

  seissol::initializers::LTS_RateAndStateFL3 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_RateAndStateFL3 *>(dynRup);
  stateVar                        = layerData.var(ConcreteLts->stateVar);
  RS_sl0                          = layerData.var(ConcreteLts->RS_sl0);
  RS_a                            = layerData.var(ConcreteLts->RS_a);

  /*a
   * Add new LTS parameter specific for this
   */
}

void seissol::dr::friction_law::VelocityWeakening::evaluate(seissol::initializers::Layer&  layerData,
                                                            seissol::initializers::DynamicRupture *dynRup,
                                                            real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                                                            real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                                                            real fullUpdateTime,
                                                            double timeWeights[CONVERGENCE_ORDER]) {
  VelocityWeakening::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
//initialize struct for in/outputs stresses
    FaultStresses faultStresses = {};

//compute stresses from Qinterpolated
    precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

/*
double RS_a;    //DISC%DynRup%RS_a  !< RS constitutive parameter "a", direct effect
double RS_sl0;  //DISC%DynRup%RS_sl0     !< Reference slip  , Dc, char. lengt scale
*/

//initialize local variables
    real LocSlip, LocSlip1, LocSlip2, LocSR1, LocSR2;
    real LocSV;
    real P_0;
    real LocP;
    real time_inc;
    real P;
    real ShTest;
    real SV0;
    real Tc;
    real coeft;
    int nSRupdates, nSVupdates;
    real SRtest;
    real tmp;
    real NR, dNR;
    real LocMu;
    real LocTracXY, LocTracXZ;

//TODO: test padded?
    for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {

      LocSlip = slip[ltsFace][iBndGP]; //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
      LocSlip1 = slipStrike[ltsFace][iBndGP]; //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
      LocSlip2 = slipDip[ltsFace][iBndGP]; //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
      LocSR1 = slipRateStrike[ltsFace][iBndGP]; //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
      LocSR2 = slipRateDip[ltsFace][iBndGP]; //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
      LocSV = stateVar[ltsFace][iBndGP];     //DISC%DynRup%StateVar(iBndGP,iFace)
      P_0 = initialStressInFaultCS[ltsFace][iBndGP][0]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];
      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {
        LocP = faultStresses.NormalStressGP[iTimeGP][iBndGP];
        time_inc = deltaT[iTimeGP];


// load traction and normal stress
        P = LocP + P_0;
        ShTest = std::sqrt(
            seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP], 2) +
            seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP], 2));

// We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
// ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

        SV0 = LocSV;    // Careful, the SV must always be corrected using SV0 and not LocSV!

// The following process is adapted from that described by Kaneko et al. (2008)
        nSRupdates = 5; //TODO: can be put outside of loop
        nSVupdates = 2;

        SlipRateMagnitude[ltsFace][iBndGP] = std::sqrt(seissol::dr::aux::power(LocSR1, 2) + seissol::dr::aux::power(LocSR2, 2)); //can be put outside of the loop

//charact. time scale Tc
        Tc = RS_sl0[ltsFace] / m_Params->rs_sr0;
// exponent
        coeft= exp(-time_inc / Tc);

        for (int j = 0; j < nSVupdates; j++) { //!This loop corrects SV values
          SlipRateMagnitude[ltsFace][iBndGP] = abs(SlipRateMagnitude[ltsFace][iBndGP]);
//exact integration assuming constant V in this loop
          LocSV= Tc * SlipRateMagnitude[ltsFace][iBndGP] * (1.0 - coeft) + coeft * SV0;
//! Newton-Raphson algorithm to determine the value of the slip rate.
//! We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
//!  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
//! In our case we equalize the values of the traction for two equations:
//!             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
//!             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
//!             where mu = mu_s + a V/(V+Vc) - b SV/(SV + Vc)
//!
          SRtest = SlipRateMagnitude[ltsFace][iBndGP];   // We use as first guess the SR value of the previous time step
          for (int i = 0; i < nSRupdates; i++) {   //!This loop corrects SR values
            tmp = m_Params->rs_f0+ RS_a[ltsFace] *SRtest/(SRtest+m_Params->rs_sr0)-m_Params->rs_b*LocSV/(LocSV+RS_sl0[ltsFace]); //=mu
            NR = - impAndEta[ltsFace].inv_eta_s * (abs(P)*tmp-ShTest)-SRtest;
            dNR          = -impAndEta[ltsFace].inv_eta_s *
                           (abs(P)*(RS_a[ltsFace]/(SRtest+m_Params->rs_sr0)-RS_a[ltsFace]*SRtest/seissol::dr::aux::power(SRtest+m_Params->rs_sr0,2))) -1.0;
            SRtest = SRtest-NR/dNR;
          }   // End nSRupdates-Loop
          tmp=0.5*(SlipRateMagnitude[ltsFace][iBndGP] + abs(SRtest));  // For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
          SlipRateMagnitude[ltsFace][iBndGP]=abs(SRtest);
        }   // End nSVupdates-Loop -  This loop corrects SV values

        LocSV    = Tc*tmp*(1-coeft) + coeft*SV0;
        tmp = 0.5 * (SlipRateMagnitude[ltsFace][iBndGP]) / m_Params->rs_sr0 * exp((m_Params->rs_f0 + m_Params->rs_b * log(m_Params->rs_sr0 * LocSV / RS_sl0[ltsFace])) / RS_a[ltsFace]);

//! Ampuero and Ben-Zion 2008 (eq. 1):
// LocMu = friction coefficient (mu_f)
// RS_f0 = static coefficient (mu_s)
// RS_a = positive coefficient, quantifying  the direct effect (alpha)
// LocSR = slip rate (V)
// RS_sr0 = characteristic velocity scale (V_c)
// RS_b = positive coefficient, quantifying  the evolution effect (beta)
// RS_sl0 = characteristic velocity scale (V_c)
        LocMu = m_Params->rs_f0 + RS_a[ltsFace] * SlipRateMagnitude[ltsFace][iBndGP] / (SlipRateMagnitude[ltsFace][iBndGP] + m_Params->rs_sr0) - m_Params->rs_b * LocSV / (LocSV + RS_sl0[ltsFace]);

// update stress change
        LocTracXY = -((initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP]) / ShTest) *LocMu * P;
        LocTracXZ = -((initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP]) / ShTest) *LocMu * P;
        LocTracXY = LocTracXY - initialStressInFaultCS[ltsFace][iBndGP][3];
        LocTracXZ = LocTracXZ - initialStressInFaultCS[ltsFace][iBndGP][5];

// Compute slip
        LocSlip = LocSlip + SlipRateMagnitude[ltsFace][iBndGP] * time_inc; // ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening

//Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
        LocSR1 = -impAndEta[ltsFace].inv_eta_s  * (LocTracXY - faultStresses.XYStressGP[iTimeGP][iBndGP]);
        LocSR2 = -impAndEta[ltsFace].inv_eta_s  * (LocTracXZ - faultStresses.XZStressGP[iTimeGP][iBndGP]);

        LocSlip1 = LocSlip1 + (LocSR1) * time_inc;
        LocSlip2 = LocSlip2 + (LocSR2) * time_inc;


//Save traction for flux computation
        faultStresses.XYTractionResultGP[iTimeGP][iBndGP] = LocTracXY;
        faultStresses.XZTractionResultGP[iTimeGP][iBndGP] = LocTracXZ;
      } //End of iTimeGP loop

      mu[ltsFace][iBndGP] = LocMu;
      slipRateStrike[ltsFace][iBndGP] = LocSR1;
      slipRateDip[ltsFace][iBndGP] = LocSR2;
      slip[ltsFace][iBndGP] = LocSlip;
      slipStrike[ltsFace][iBndGP] = LocSlip1;
      slipDip[ltsFace][iBndGP] = LocSlip2;
      stateVar[ltsFace][iBndGP] = LocSV;
      tractionXY[ltsFace][iBndGP] = LocTracXY;
      tractionXZ[ltsFace][iBndGP] = LocTracXZ;
    }//End of iBndGP-loop

// output rupture front
// outside of iTimeGP loop in order to safe an 'if' in a loop
// this way, no subtimestep resolution possible
    saveRuptureFrontOutput(ltsFace);

//output peak slip rate
    savePeakSlipRateOutput(ltsFace);

//save stresses in imposedState
    postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, timeWeights, ltsFace);

  }//End of Loop over Faces

}//End of Function evaluate
