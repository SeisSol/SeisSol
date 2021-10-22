#include "ImposedSlipRates.h"

namespace seissol::dr::friction_law {
void ImposedSlipRates::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                          seissol::initializers::DynamicRupture* dynRup,
                                          real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

  seissol::initializers::LTS_ImposedSlipRatesFL33* ConcreteLts =
      dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesFL33*>(dynRup);
  nucleationStressInFaultCS = layerData.var(ConcreteLts->nucleationStressInFaultCS);
  averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
}

void ImposedSlipRates::evaluate(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real fullUpdateTime,
    double timeWeights[CONVERGENCE_ORDER]) {

  copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
    // initialize struct for in/outputs stresses
    FaultStresses faultStresses = {};

    // declare local variables
    std::array<real, numOfPointsPadded> tmpSlip{0};
    real tn = fullUpdateTime;
    real time_inc;
    real Gnuc = 0.0;

    // compute stresses from Qinterpolated
    precomputeStressFromQInterpolated(
        faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

    for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) { // loop over time steps
      time_inc = deltaT[iTimeGP];
      tn = tn + time_inc;
      Gnuc = Calc_SmoothStepIncrement(tn, time_inc) / time_inc;

      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
        //! EQN%NucleationStressInFaultCS (1 and 2) contains the slip in FaultCS
        faultStresses.XYTractionResultGP[iTimeGP][iBndGP] =
            faultStresses.XYStressGP[iTimeGP][iBndGP] -
            impAndEta[ltsFace].eta_s * nucleationStressInFaultCS[ltsFace][iBndGP][0] * Gnuc;
        faultStresses.XZTractionResultGP[iTimeGP][iBndGP] =
            faultStresses.XZStressGP[iTimeGP][iBndGP] -
            impAndEta[ltsFace].eta_s * nucleationStressInFaultCS[ltsFace][iBndGP][1] * Gnuc;
        slipRateStrike[ltsFace][iBndGP] = nucleationStressInFaultCS[ltsFace][iBndGP][0] * Gnuc;
        slipRateDip[ltsFace][iBndGP] = nucleationStressInFaultCS[ltsFace][iBndGP][1] * Gnuc;
        SlipRateMagnitude[ltsFace][iBndGP] =
            std::sqrt(seissol::dr::aux::power(slipRateStrike[ltsFace][iBndGP], 2) +
                      seissol::dr::aux::power(slipRateDip[ltsFace][iBndGP], 2));

        //! Update slip
        slipStrike[ltsFace][iBndGP] += slipRateStrike[ltsFace][iBndGP] * time_inc;
        slipDip[ltsFace][iBndGP] += slipRateDip[ltsFace][iBndGP] * time_inc;
        slip[ltsFace][iBndGP] += SlipRateMagnitude[ltsFace][iBndGP] * time_inc;
        tmpSlip[iBndGP] += SlipRateMagnitude[ltsFace][iBndGP] * time_inc;

        tractionXY[ltsFace][iBndGP] = faultStresses.XYTractionResultGP[iTimeGP][iBndGP];
        tractionXZ[ltsFace][iBndGP] = faultStresses.XYTractionResultGP[iTimeGP][iBndGP];
      }
    }
    // output rupture front
    // outside of iTimeGP loop in order to safe an 'if' in a loop
    // this way, no subtimestep resolution possible
    saveRuptureFrontOutput(ltsFace);

    // output peak slip rate
    savePeakSlipRateOutput(ltsFace);

    //---compute and store slip to determine the magnitude of an earthquake ---
    //    to this end, here the slip is computed and averaged per element
    //    in calc_seissol.f90 this value will be multiplied by the element surface
    //    and an output happened once at the end of the simulation
    saveAverageSlipOutput(tmpSlip, ltsFace);

    // save stresses in imposedState
    postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace],
                                         QInterpolatedMinus[ltsFace],
                                         faultStresses,
                                         timeWeights,
                                         ltsFace);
  } // End of Loop over Faces
}
} // namespace seissol::dr::friction_law
