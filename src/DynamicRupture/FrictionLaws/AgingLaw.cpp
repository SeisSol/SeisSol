#include "AgingLaw.h"
namespace seissol::dr::friction_law {
real AgingLaw::calcStateVariableHook(real stateVariable, real tmp, real time_inc, real RS_sl0) {
  return stateVariable * std::exp(-tmp * time_inc / RS_sl0) +
         RS_sl0 / tmp * (1.0 - std::exp(-tmp * time_inc / RS_sl0));
}

void AgingLaw::evaluate(
    initializers::Layer& layerData,
    initializers::DynamicRupture* dynRup,
    real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real fullUpdateTime,
    double timeWeights[CONVERGENCE_ORDER]) {
  initializers::LTS_RateAndStateFL3* ConcreteLts =
      dynamic_cast<initializers::LTS_RateAndStateFL3*>(dynRup);

  model::IsotropicWaveSpeeds* waveSpeedsPlus = layerData.var(ConcreteLts->waveSpeedsPlus);
  model::IsotropicWaveSpeeds* waveSpeedsMinus = layerData.var(ConcreteLts->waveSpeedsMinus);
  real(*initialStressInFaultCS)[numPaddedPoints][6] =
      layerData.var(ConcreteLts->initialStressInFaultCS);
  real(*cohesion)[numPaddedPoints] = layerData.var(ConcreteLts->cohesion);

  real* RS_a = layerData.var(ConcreteLts->RS_a);
  real* RS_sl0 = layerData.var(ConcreteLts->RS_sl0);
  real* RS_sr0 = layerData.var(ConcreteLts->RS_sr0);

  real(*mu)[numPaddedPoints] = layerData.var(ConcreteLts->mu);
  real(*slip)[numPaddedPoints] = layerData.var(ConcreteLts->slip);
  real(*slip1)[numPaddedPoints] = layerData.var(ConcreteLts->slipStrike);
  real(*slip2)[numPaddedPoints] = layerData.var(ConcreteLts->slipDip);
  real(*slipRate1)[numPaddedPoints] = layerData.var(ConcreteLts->slipRateStrike);
  real(*slipRate2)[numPaddedPoints] = layerData.var(ConcreteLts->slipRateDip);
  real(*StateVar)[numPaddedPoints] = layerData.var(ConcreteLts->stateVar);

  real(*tracXY)[numPaddedPoints] = layerData.var(ConcreteLts->tractionXY);
  real(*tracXZ)[numPaddedPoints] = layerData.var(ConcreteLts->tractionXZ);

  // loop parameter are fixed, not variable??
  constexpr unsigned int nSRupdates{5}, nSVupdates{2};

#ifdef _OPENMP
#pragma omp parallel for schedule(static) // private(QInterpolatedPlus,QInterpolatedMinus)
#endif
  for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {

    FaultStresses faultStresses = {};

    precomputeStressFromQInterpolated(
        faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

      // Find variables at given fault node
      real localSlip = slip[ltsFace][pointIndex];              // Slip path
      real localSlip1 = slip1[ltsFace][pointIndex];            // Slip along direction 1
      real localSlip2 = slip2[ltsFace][pointIndex];            // Slip along direction 2
      real localSR1 = slipRate1[ltsFace][pointIndex];          // Slip rate along direction 1
      real localSR2 = slipRate2[ltsFace][pointIndex];          // Slip rate along direction 2
      real localStateVariable = StateVar[ltsFace][pointIndex]; // State Variable
      real localCohesion = cohesion[ltsFace][pointIndex];      // cohesion: should be negative since
                                                          // negative normal stress is compression
      real initialPressure = initialStressInFaultCS[ltsFace][pointIndex][0]; // initial pressure

      real localMu = 0;
      real localTractionXY = 0;
      real localTractionXZ = 0;

      for (int timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
        real localPressure = faultStresses.NormalStressGP[timeIndex][pointIndex];
        real timeIncrement = deltaT[timeIndex];

        // load traction and normal stress
        real pressure = localPressure + initialPressure;

        real stressXY = initialStressInFaultCS[ltsFace][pointIndex][3] +
                        faultStresses.XYStressGP[timeIndex][pointIndex];
        real stressXZ = initialStressInFaultCS[ltsFace][pointIndex][5] +
                        faultStresses.XZStressGP[timeIndex][pointIndex];
        real totalShearStressYZ = std::sqrt(std::pow(stressXY, 2) + std::pow(stressXZ, 2));

        // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
        // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

        // Careful, the stateVariable must always be corrected using stateVariable and not
        // localStateVariable!
        real stateVariable = localStateVariable;

        // The following process is adapted from that described by Kaneko et al. (2008)
        slipRateMagnitude[ltsFace][pointIndex] =
            std::sqrt(std::pow(localSR1, 2) + std::pow(localSR2, 2));
        real tmp = std::fabs(slipRateMagnitude[ltsFace][pointIndex]);

        // This loop corrects StateVariable values
        for (unsigned int j = 0; j < nSVupdates; j++) {
          slipRateMagnitude[ltsFace][pointIndex] =
              std::fabs(slipRateMagnitude[ltsFace][pointIndex]);

          // FL= 3 aging law and FL=4 slip law
          localStateVariable =
              calcStateVariableHook(stateVariable, tmp, timeIncrement, RS_sl0[ltsFace]);

          // Newton-Raphson algorithm to determine the value of the slip rate.
          // We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g ,
          // which has
          //  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i
          //  / dNR_i ).
          // In our case we equalize the values of the traction for two equations:
          //             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
          //             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
          //               where mu=a*asinh(SR/2/SR0*exp((F0+b*log(SR0*SV/L))/a (eq. 2a of Lapusta
          //               and Rice (2003))

          // SRtest: We use as first guess the  SR value of the previous time step
          real slipRateGuess = slipRateMagnitude[ltsFace][pointIndex];

          for (unsigned int i = 0; i < nSRupdates; i++) { // This loop corrects SR values
            tmp = 0.5 / RS_sr0[ltsFace] *
                  std::exp((m_Params->rs_f0 +
                            m_Params->rs_b *
                                std::log(RS_sr0[ltsFace] * localStateVariable / RS_sl0[ltsFace])) /
                           RS_a[ltsFace]);
            real tmp2 = tmp * slipRateGuess;
            // TODO: author before me: not sure if ShTest=TotalShearStressYZ should be + or -...
            real NR = -(1.0 / waveSpeedsPlus->sWaveVelocity / waveSpeedsPlus->density +
                        1.0 / waveSpeedsMinus->sWaveVelocity / waveSpeedsMinus->density) *
                          (std::fabs(pressure) * RS_a[ltsFace] *
                               std::log(tmp2 + std::sqrt(std::pow(tmp2, 2) + 1.0)) -
                           totalShearStressYZ) -
                      slipRateGuess;

            real dNR =
                -(1.0 / waveSpeedsPlus->sWaveVelocity / waveSpeedsPlus->density +
                  1.0 / waveSpeedsMinus->sWaveVelocity / waveSpeedsMinus->density) *
                    (std::fabs(pressure) * RS_a[ltsFace] / std::sqrt(1 + std::pow(tmp2, 2)) * tmp) -
                1.0;
            // no ABS needed around NR/dNR at least for aging law
            real slipRateGuess = std::fabs(slipRateGuess - NR / dNR);
          }

          // For the next SV update, use the mean slip rate between the initial guess and the one
          // found (Kaneko 2008, step 6)
          tmp = 0.5 * (slipRateMagnitude[ltsFace][pointIndex] + std::fabs(slipRateGuess));

          slipRateMagnitude[ltsFace][pointIndex] = std::fabs(slipRateGuess);
        } // End SV-Loop

        // FL= 3 aging law and FL=4 slip law
        localStateVariable =
            calcStateVariableHook(stateVariable, tmp, timeIncrement, RS_sl0[ltsFace]);

        // TODO: reused calc from above -> simplify
        tmp = 0.5 * (slipRateMagnitude[ltsFace][pointIndex]) / RS_sr0[ltsFace] *
              std::exp((m_Params->rs_f0 +
                        m_Params->rs_b *
                            std::log(RS_sr0[ltsFace] * localStateVariable / RS_sl0[ltsFace])) /
                       RS_a[ltsFace]);

        localMu = RS_a[ltsFace] * std::log(tmp + std::sqrt(std::pow(tmp, 2) + 1.0));

        // 2D:
        // LocTrac  = -(ABS(S_0)-LocMu*(LocP+P_0))*(S_0/ABS(S_0))
        // LocTrac  = ABS(LocTrac)*(-SignSR)  !!! line commented as it leads NOT to correct results
        // update stress change
        localTractionXY = -((initialStressInFaultCS[ltsFace][pointIndex][3] +
                             faultStresses.XYStressGP[pointIndex][timeIndex]) /
                            totalShearStressYZ) *
                          (localMu * pressure + std::fabs(localCohesion));
        localTractionXZ = -((initialStressInFaultCS[ltsFace][pointIndex][5] +
                             faultStresses.XZStressGP[pointIndex][timeIndex]) /
                            totalShearStressYZ) *
                          (localMu * pressure + std::fabs(localCohesion));
        localTractionXY = localTractionXY - initialStressInFaultCS[ltsFace][pointIndex][3];
        localTractionXZ = localTractionXZ - initialStressInFaultCS[ltsFace][pointIndex][5];

        // Compute slip
        // ABS of LocSR removed as it would be the accumulated slip that is usually not needed in
        // the solver, see linear slip weakening
        localSlip = localSlip + (slipRateMagnitude[ltsFace][pointIndex]) * timeIncrement;

        // Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused
        // by a free surface!)
        localSR1 = -(1.0 / (waveSpeedsPlus->sWaveVelocity * waveSpeedsPlus->density) +
                     1.0 / (waveSpeedsMinus->sWaveVelocity * waveSpeedsMinus->density)) *
                   (localTractionXY - faultStresses.XYStressGP[timeIndex][pointIndex]);
        localSR2 = -(1.0 / (waveSpeedsPlus->sWaveVelocity * waveSpeedsPlus->density) +
                     1.0 / (waveSpeedsMinus->sWaveVelocity * waveSpeedsMinus->density)) *
                   (localTractionXZ - faultStresses.XZStressGP[timeIndex][pointIndex]);

        localSlip1 = localSlip1 + localSR1 * timeIncrement;
        localSlip2 = localSlip2 + localSR2 * timeIncrement;

        // Save traction for flux computation
        faultStresses.XYTractionResultGP[timeIndex][pointIndex] = localTractionXY;
        faultStresses.XZTractionResultGP[timeIndex][pointIndex] = localTractionXZ;
      } // End of timeIndex- loop

      mu[ltsFace][pointIndex] = localMu;
      slipRate1[ltsFace][pointIndex] = localSR1;
      slipRate2[ltsFace][pointIndex] = localSR2;
      slip[ltsFace][pointIndex] = localSlip;
      slip1[ltsFace][pointIndex] = localSlip1;
      slip2[ltsFace][pointIndex] = localSlip2;
      StateVar[ltsFace][pointIndex] = localStateVariable;
      tracXY[ltsFace][pointIndex] = localTractionXY;
      tracXZ[ltsFace][pointIndex] = localTractionXZ;

    } // End of pointIndex-loop

    // output rupture front
    // outside of timeIndex loop in order to safe an 'if' in a loop
    // this way, no subtimestep resolution possible
    saveRuptureFrontOutput(ltsFace);

    savePeakSlipRateOutput(ltsFace);

    postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace],
                                         QInterpolatedMinus[ltsFace],
                                         faultStresses,
                                         timeWeights,
                                         ltsFace);
  } // end face-loop
} // end evaluate function

} // namespace seissol::dr::friction_law