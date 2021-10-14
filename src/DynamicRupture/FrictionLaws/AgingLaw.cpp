#include "AgingLaw.h"

real seissol::dr::friction_law::AgingLaw::calcStateVariableHook(real SV0, real tmp, real time_inc, real RS_sl0) {
  return SV0*exp(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-exp(-tmp*time_inc/RS_sl0));
}

void seissol::dr::friction_law::AgingLaw::evaluate(seissol::initializers::Layer&  layerData,
                                                   seissol::initializers::DynamicRupture *dynRup,
                                                   real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                                                   real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                                                   real fullUpdateTime,
                                                   real timeWeights[CONVERGENCE_ORDER]) {
  seissol::initializers::LTS_RateAndStateFL3 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_RateAndStateFL3 *>(dynRup);

  seissol::model::IsotropicWaveSpeeds *waveSpeedsPlus                           = layerData.var(ConcreteLts->waveSpeedsPlus);
  seissol::model::IsotropicWaveSpeeds *waveSpeedsMinus                          = layerData.var(ConcreteLts->waveSpeedsMinus);
  real                    (*initialStressInFaultCS)[numOfPointsPadded][6]       = layerData.var(ConcreteLts->initialStressInFaultCS);
  real                    (*cohesion)[numOfPointsPadded]                        = layerData.var(ConcreteLts->cohesion);

  real*                   RS_a                                                  = layerData.var(ConcreteLts->RS_a);
  real*                   RS_sl0                                                = layerData.var(ConcreteLts->RS_sl0);
  real*                   RS_sr0                                                = layerData.var(ConcreteLts->RS_sr0);

  real                    (*mu)[numOfPointsPadded]                              = layerData.var(ConcreteLts->mu);
  real                    (*slip)[numOfPointsPadded]                            = layerData.var(ConcreteLts->slip);
  real                    (*slip1)[numOfPointsPadded]                           = layerData.var(ConcreteLts->slipStrike);
  real                    (*slip2)[numOfPointsPadded]                           = layerData.var(ConcreteLts->slipDip);
  real                    (*slipRate1)[numOfPointsPadded]                       = layerData.var(ConcreteLts->slipRateStrike);
  real                    (*slipRate2)[numOfPointsPadded]                       = layerData.var(ConcreteLts->slipRateDip);
  real                    (*StateVar)[numOfPointsPadded]                        = layerData.var(ConcreteLts->stateVar);

  real                    (*tracXY)[numOfPointsPadded]                          = layerData.var(ConcreteLts->tractionXY);
  real                    (*tracXZ)[numOfPointsPadded]                          = layerData.var(ConcreteLts->tractionXZ);


//loop parameter are fixed, not variable??
  unsigned int nSRupdates, nSVupdates;
  nSRupdates = 5;
  nSVupdates = 2;


#ifdef _OPENMP
#pragma omp parallel for schedule(static) //private(QInterpolatedPlus,QInterpolatedMinus)
#endif
  for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {

    FaultStresses faultStresses = {};

    precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

    real LocSlip, LocSlip1, LocSlip2, LocSR1, LocSR2, LocSV, LocCohesion, P_0, LocP, time_inc, P, TotalShearStressYZ, SV0, tmp,tmp2, SlipRateGuess, NR, dNR, LocMu;
    real LocTracXY, LocTracXZ;

    for(int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      LocSlip = slip[ltsFace][iBndGP]; //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
      LocSlip1 = slip1[ltsFace][iBndGP]; //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
      LocSlip2 = slip2[ltsFace][iBndGP]; //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
      LocSR1 = slipRate1[ltsFace][iBndGP]; //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
      LocSR2 = slipRate2[ltsFace][iBndGP]; //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
      LocSV = StateVar[ltsFace][iBndGP];     //DISC%DynRup%StateVar(iBndGP,iFace)
      LocCohesion = cohesion[ltsFace][iBndGP]; //DISC%DynRup%cohesion(iBndGP,iFace)          // !< cohesion at given fault node  (should be negative since negative normal stress is compression)
      P_0 = initialStressInFaultCS[ltsFace][iBndGP][0]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {
        LocP = faultStresses.NormalStressGP[iTimeGP][iBndGP];
        time_inc = deltaT[iTimeGP];

//SignSR1   = SIGN(1.0,LocSR1)                    ! Gets the sign of the slip rate
//SignSR2   = SIGN(1.0,LocSR2)                    ! Gets the sign of the slip rate

// load traction and normal stress
        P = LocP + P_0;

        TotalShearStressYZ = std::sqrt(
            seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP], 2) +
            seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP], 2));

// We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
// ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

        SV0 = LocSV;    // Careful, the SV must always be corrected using SV0 and not LocSV!

// The following process is adapted from that described by Kaneko et al. (2008)
        SlipRateMagnitude[ltsFace][iBndGP]      = std::sqrt(seissol::dr::aux::power(LocSR1, 2) + seissol::dr::aux::power(LocSR2, 2));
        tmp        = fabs(SlipRateMagnitude[ltsFace][iBndGP]);

        for(unsigned int j = 0; j < nSVupdates; j++){ //!This loop corrects SV values
          SlipRateMagnitude[ltsFace][iBndGP] = fabs(SlipRateMagnitude[ltsFace][iBndGP]);

// FL= 3 aging law and FL=4 slip law
          LocSV = calcStateVariableHook( SV0,  tmp,  time_inc,  RS_sl0[ltsFace]);

// Newton-Raphson algorithm to determine the value of the slip rate.
// We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
//  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
// In our case we equalize the values of the traction for two equations:
//             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
//             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
//               where mu=a*asinh(SR/2/SR0*exp((F0+b*log(SR0*SV/L))/a (eq. 2a of Lapusta and Rice (2003))

          SlipRateGuess = SlipRateMagnitude[ltsFace][iBndGP];   // SRtest: We use as first guess the SR value of the previous time step

          for(unsigned int i = 0; i < nSRupdates; i++){   //!This loop corrects SR values
            tmp          = 0.5 / RS_sr0[ltsFace] * exp((m_Params->rs_f0 + m_Params->rs_b * log(RS_sr0[ltsFace] * LocSV / RS_sl0[ltsFace]) ) / RS_a[ltsFace]);
            tmp2         = tmp * SlipRateGuess;
            NR           = -(1.0/waveSpeedsPlus->sWaveVelocity/waveSpeedsPlus->density+1.0/waveSpeedsMinus->sWaveVelocity/waveSpeedsMinus->density) *
                           (fabs(P) * RS_a[ltsFace] * log(tmp2 + sqrt(seissol::dr::aux::power(tmp2, 2) + 1.0)) - TotalShearStressYZ) - SlipRateGuess;    //!TODO: author before me: not sure if ShTest=TotalShearStressYZ should be + or -...
            dNR          = -(1.0/waveSpeedsPlus->sWaveVelocity/waveSpeedsPlus->density+1.0/waveSpeedsMinus->sWaveVelocity/waveSpeedsMinus->density) *
                           (fabs(P) * RS_a[ltsFace] / sqrt(1 + seissol::dr::aux::power(tmp2, 2)) * tmp) - 1.0;
            SlipRateGuess = fabs(SlipRateGuess-NR/dNR);             // no ABS needed around NR/dNR at least for aging law
          }   // End
          tmp=0.5*(SlipRateMagnitude[ltsFace][iBndGP] + fabs(SlipRateGuess));  //! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
          SlipRateMagnitude[ltsFace][iBndGP]=fabs(SlipRateGuess);
        }   // End SV-Loop

// FL= 3 aging law and FL=4 slip law
        LocSV= calcStateVariableHook( SV0,  tmp,  time_inc,  RS_sl0[ltsFace]);

//TODO: reused calc from above -> simplify
        tmp  = 0.5 * ( SlipRateMagnitude[ltsFace][iBndGP]) / RS_sr0[ltsFace] * exp((m_Params->rs_f0 + m_Params->rs_b * log(RS_sr0[ltsFace] * LocSV / RS_sl0[ltsFace])) / RS_a[ltsFace]);

        LocMu    = RS_a[ltsFace] * log(tmp + sqrt(seissol::dr::aux::power(tmp, 2) + 1.0));

// 2D:
// LocTrac  = -(ABS(S_0)-LocMu*(LocP+P_0))*(S_0/ABS(S_0))
// LocTrac  = ABS(LocTrac)*(-SignSR)  !!! line commented as it leads NOT to correct results
// update stress change
        LocTracXY = -((initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iBndGP][iTimeGP]) / TotalShearStressYZ) * (LocMu * P + fabs(LocCohesion));
        LocTracXZ = -((initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iBndGP][iTimeGP]) / TotalShearStressYZ) * (LocMu * P + fabs(LocCohesion));
        LocTracXY = LocTracXY - initialStressInFaultCS[ltsFace][iBndGP][3];
        LocTracXZ = LocTracXZ - initialStressInFaultCS[ltsFace][iBndGP][5];

// Compute slip
        LocSlip   = LocSlip  + ( SlipRateMagnitude[ltsFace][iBndGP]) * time_inc; // ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening

//Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
        LocSR1     = -(1.0/(waveSpeedsPlus->sWaveVelocity*waveSpeedsPlus->density)+1.0/(waveSpeedsMinus->sWaveVelocity*waveSpeedsMinus->density))*(LocTracXY-faultStresses.XYStressGP[iTimeGP][iBndGP]);
        LocSR2     = -(1.0/(waveSpeedsPlus->sWaveVelocity*waveSpeedsPlus->density)+1.0/(waveSpeedsMinus->sWaveVelocity*waveSpeedsMinus->density))*(LocTracXZ-faultStresses.XZStressGP[iTimeGP][iBndGP]);

        LocSlip1   = LocSlip1  + (LocSR1)*time_inc;
        LocSlip2   = LocSlip2  + (LocSR2)*time_inc;

//LocSR1     = SignSR1*ABS(LocSR1)
//LocSR2     = SignSR2*ABS(LocSR2)

//Save traction for flux computation
        faultStresses.XYTractionResultGP[iTimeGP][iBndGP] = LocTracXY;
        faultStresses.XZTractionResultGP[iTimeGP][iBndGP] = LocTracXZ;
      }//End of iTimeGP- loop

      mu[ltsFace][iBndGP]        = LocMu;
      slipRate1[ltsFace][iBndGP] = LocSR1;
      slipRate2[ltsFace][iBndGP] = LocSR2;
      slip[ltsFace][iBndGP]      = LocSlip;
      slip1[ltsFace][iBndGP]     = LocSlip1;
      slip2[ltsFace][iBndGP]     = LocSlip2;
      StateVar[ltsFace][iBndGP]  = LocSV;
      tracXY[ltsFace][iBndGP]    = LocTracXY;
      tracXZ[ltsFace][iBndGP]    = LocTracXZ;

    }//End of iBndGP-loop

// output rupture front
// outside of iTimeGP loop in order to safe an 'if' in a loop
// this way, no subtimestep resolution possible
    saveRuptureFrontOutput(ltsFace);

    savePeakSlipRateOutput(ltsFace);


    postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace],
                                         faultStresses, timeWeights, ltsFace);
  } //end face-loop
} //end evaluate function

