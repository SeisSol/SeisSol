#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
void BaseFrictionLaw::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                         seissol::initializers::DynamicRupture* dynRup,
                                         real fullUpdateTime) {
  impAndEta = layerData.var(dynRup->impAndEta);
  initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS);
  mu = layerData.var(dynRup->mu);
  slip = layerData.var(dynRup->slip);
  slipStrike = layerData.var(dynRup->slipStrike);
  slipDip = layerData.var(dynRup->slipDip);
  slipRateMagnitude = layerData.var(dynRup->slipRateMagnitude);
  slipRateStrike = layerData.var(dynRup->slipRateStrike);
  slipRateDip = layerData.var(dynRup->slipRateDip);
  ruptureTime = layerData.var(dynRup->ruptureTime);
  ruptureFront = layerData.var(dynRup->ruptureFront);
  peakSlipRate = layerData.var(dynRup->peakSlipRate);
  tractionXY = layerData.var(dynRup->tractionXY);
  tractionXZ = layerData.var(dynRup->tractionXZ);
  imposedStatePlus = layerData.var(dynRup->imposedStatePlus);
  imposedStateMinus = layerData.var(dynRup->imposedStateMinus);
  m_fullUpdateTime = fullUpdateTime;
}

void BaseFrictionLaw::precomputeStressFromQInterpolated(
    FaultStresses& faultStresses,
    real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    unsigned int ltsFace) {
  // this initialization of the kernel could be moved to the initializer,
  // since all inputs outside the j-loop are time independent
  // set inputParam could be extendent for this
  // the kernel then could be a class attribute (but be careful of race conditions since this is
  // computed in parallel!!)
  dynamicRupture::kernel::StressFromQInterpolated StressFromQInterpolatedKrnl;
  StressFromQInterpolatedKrnl.eta_p = impAndEta[ltsFace].eta_p;
  StressFromQInterpolatedKrnl.eta_s = impAndEta[ltsFace].eta_s;
  StressFromQInterpolatedKrnl.inv_Zp = impAndEta[ltsFace].inv_Zp;
  StressFromQInterpolatedKrnl.inv_Zs = impAndEta[ltsFace].inv_Zs;
  StressFromQInterpolatedKrnl.inv_Zp_neig = impAndEta[ltsFace].inv_Zp_neig;
  StressFromQInterpolatedKrnl.inv_Zs_neig = impAndEta[ltsFace].inv_Zs_neig;
  StressFromQInterpolatedKrnl.select0 = init::select0::Values;
  StressFromQInterpolatedKrnl.select3 = init::select3::Values;
  StressFromQInterpolatedKrnl.select5 = init::select5::Values;
  StressFromQInterpolatedKrnl.select6 = init::select6::Values;
  StressFromQInterpolatedKrnl.select7 = init::select7::Values;
  StressFromQInterpolatedKrnl.select8 = init::select8::Values;

  for (int j = 0; j < CONVERGENCE_ORDER; j++) {
    StressFromQInterpolatedKrnl.QInterpolatedMinus = QInterpolatedMinus[j];
    StressFromQInterpolatedKrnl.QInterpolatedPlus = QInterpolatedPlus[j];
    StressFromQInterpolatedKrnl.NorStressGP = faultStresses.NormalStressGP[j];
    StressFromQInterpolatedKrnl.XYStressGP = faultStresses.XYStressGP[j];
    StressFromQInterpolatedKrnl.XZStressGP = faultStresses.XZStressGP[j];
    // Carsten Uphoff Thesis: EQ.: 4.53
    StressFromQInterpolatedKrnl.execute();
  }

  static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                "Different number of quadrature points?");
}

void BaseFrictionLaw::postcomputeImposedStateFromNewStress(
    real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const FaultStresses& faultStresses,
    double timeWeights[CONVERGENCE_ORDER],
    unsigned int ltsFace) {
  // this initialization of the kernel could be moved to the initializer
  // set inputParam could be extendent for this (or create own function)
  // the kernel then could be a class attribute and following values are only set once
  //(but be careful of race conditions since this is computed in parallel for each face!!)
  dynamicRupture::kernel::ImposedStateFromNewStress ImposedStateFromNewStressKrnl;
  ImposedStateFromNewStressKrnl.select0 = init::select0::Values;
  ImposedStateFromNewStressKrnl.select3 = init::select3::Values;
  ImposedStateFromNewStressKrnl.select5 = init::select5::Values;
  ImposedStateFromNewStressKrnl.select6 = init::select6::Values;
  ImposedStateFromNewStressKrnl.select7 = init::select7::Values;
  ImposedStateFromNewStressKrnl.select8 = init::select8::Values;
  ImposedStateFromNewStressKrnl.inv_Zs = impAndEta[ltsFace].inv_Zs;
  ImposedStateFromNewStressKrnl.inv_Zs_neig = impAndEta[ltsFace].inv_Zs_neig;
  ImposedStateFromNewStressKrnl.inv_Zp = impAndEta[ltsFace].inv_Zp;
  ImposedStateFromNewStressKrnl.inv_Zp_neig = impAndEta[ltsFace].inv_Zp_neig;

  // set imposed state to zero
  for (unsigned int i = 0; i < tensor::QInterpolated::size(); i++) {
    imposedStatePlus[ltsFace][i] = 0;
    imposedStateMinus[ltsFace][i] = 0;
  }
  ImposedStateFromNewStressKrnl.imposedStatePlus = imposedStatePlus[ltsFace];
  ImposedStateFromNewStressKrnl.imposedStateMinus = imposedStateMinus[ltsFace];

  for (int j = 0; j < CONVERGENCE_ORDER; j++) {
    ImposedStateFromNewStressKrnl.NorStressGP = faultStresses.NormalStressGP[j];
    ImposedStateFromNewStressKrnl.TractionGP_XY = faultStresses.XYTractionResultGP[j];
    ImposedStateFromNewStressKrnl.TractionGP_XZ = faultStresses.XZTractionResultGP[j];
    ImposedStateFromNewStressKrnl.timeWeights = timeWeights[j];
    ImposedStateFromNewStressKrnl.QInterpolatedMinus = QInterpolatedMinus[j];
    ImposedStateFromNewStressKrnl.QInterpolatedPlus = QInterpolatedPlus[j];
    // Carsten Uphoff Thesis: EQ.: 4.60
    ImposedStateFromNewStressKrnl.execute();
  }
}

/*
 * https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
real BaseFrictionLaw::calcSmoothStepIncrement(real currentTime, real dt) {
  real gNuc;
  real prevTime;
  if (currentTime > 0.0 && currentTime <= drParameters.t_0) {
    gNuc = calcSmoothStep(currentTime);
    prevTime = currentTime - dt;
    if (prevTime > 0.0) {
      gNuc = gNuc - calcSmoothStep(prevTime);
    }
  } else {
    gNuc = 0.0;
  }
  return gNuc;
}

/*
 * https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
real BaseFrictionLaw::calcSmoothStep(real currentTime) {
  real gNuc;
  if (currentTime <= 0) {
    gNuc = 0.0;
  } else {
    if (currentTime < drParameters.t_0) {
      gNuc = std::exp(std::pow(currentTime - drParameters.t_0, 2) /
                      (currentTime * (currentTime - 2.0 * drParameters.t_0)));
    } else {
      gNuc = 1.0;
    }
  }
  return gNuc;
}

void BaseFrictionLaw::saveRuptureFrontOutput(unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    constexpr real ruptureFrontThreshold = 0.001;
    if (ruptureFront[ltsFace][pointIndex] &&
        slipRateMagnitude[ltsFace][pointIndex] > ruptureFrontThreshold) {
      ruptureTime[ltsFace][pointIndex] = m_fullUpdateTime;
      ruptureFront[ltsFace][pointIndex] = false;
    }
  }
}

void BaseFrictionLaw::savePeakSlipRateOutput(unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    if (slipRateMagnitude[ltsFace][pointIndex] > peakSlipRate[ltsFace][pointIndex]) {
      peakSlipRate[ltsFace][pointIndex] = slipRateMagnitude[ltsFace][pointIndex];
    }
  }
}

void BaseFrictionLaw::saveAverageSlipOutput(std::array<real, numPaddedPoints>& tmpSlip,
                                            unsigned int ltsFace) {
  real sum_tmpSlip = 0;
  if (drParameters.isMagnitudeOutputOn) {
    for (int pointIndex = 0; pointIndex < numberOfPoints; pointIndex++)
      sum_tmpSlip += tmpSlip[pointIndex];
    averagedSlip[ltsFace] = averagedSlip[ltsFace] + sum_tmpSlip / numberOfPoints;
  }
}

void BaseFrictionLaw::computeDeltaT(double* timePoints) {
  deltaT[0] = timePoints[0];
  for (int timeIndex = 1; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }
  // to fill last segment of Gaussian integration
  deltaT[CONVERGENCE_ORDER - 1] = deltaT[CONVERGENCE_ORDER - 1] + deltaT[0];
}
} // namespace seissol::dr::friction_law
