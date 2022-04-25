#include "FrictionSolver.h"

namespace seissol::dr::friction_law {

void FrictionSolver::computeDeltaT(double timePoints[CONVERGENCE_ORDER]) {
  deltaT[0] = timePoints[0];
  for (unsigned timeIndex = 1; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }
  // to fill last segment of Gaussian integration
  deltaT[CONVERGENCE_ORDER - 1] = deltaT[CONVERGENCE_ORDER - 1] + deltaT[0];
}

void FrictionSolver::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                        seissol::initializers::DynamicRupture* dynRup,
                                        real fullUpdateTime) {
  impAndEta = layerData.var(dynRup->impAndEta);
  initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS);
  mu = layerData.var(dynRup->mu);
  accumulatedSlipMagnitude = layerData.var(dynRup->accumulatedSlipMagnitude);
  slip1 = layerData.var(dynRup->slip1);
  slip2 = layerData.var(dynRup->slip2);
  slipRateMagnitude = layerData.var(dynRup->slipRateMagnitude);
  slipRate1 = layerData.var(dynRup->slipRate1);
  slipRate2 = layerData.var(dynRup->slipRate2);
  ruptureTime = layerData.var(dynRup->ruptureTime);
  ruptureTimePending = layerData.var(dynRup->ruptureTimePending);
  peakSlipRate = layerData.var(dynRup->peakSlipRate);
  traction1 = layerData.var(dynRup->traction1);
  traction2 = layerData.var(dynRup->traction2);
  imposedStatePlus = layerData.var(dynRup->imposedStatePlus);
  imposedStateMinus = layerData.var(dynRup->imposedStateMinus);
  mFullUpdateTime = fullUpdateTime;
  averagedSlip = layerData.var(dynRup->averagedSlip);
  dynStressTime = layerData.var(dynRup->dynStressTime);
  dynStressTimePending = layerData.var(dynRup->dynStressTimePending);
  qInterpolatedPlus = layerData.var(dynRup->qInterpolatedPlus);
  qInterpolatedMinus = layerData.var(dynRup->qInterpolatedMinus);
}

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp declare target
#endif // ACL_DEVICE_OFFLOAD
void FrictionSolver::precomputeStressFromQInterpolated(FaultStresses& faultStresses,
                                                       unsigned int ltsFace) {

  static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                "Different number of quadrature points?");

  // this initialization of the kernel could be moved to the initializer,
  // since all inputs outside the j-loop are time independent
  // set inputParam could be extendent for this
  // the kernel then could be a class attribute (but be careful of race conditions since this is
  // computed in parallel!!)

  auto etaP = impAndEta[ltsFace].etaP;
  auto etaS = impAndEta[ltsFace].etaS;
  auto invZp = impAndEta[ltsFace].invZp;
  auto invZs = impAndEta[ltsFace].invZs;
  auto invZpNeig = impAndEta[ltsFace].invZpNeig;
  auto invZsNeig = impAndEta[ltsFace].invZsNeig;

  using QInterpolatedShapeT =
      real(*)[CONVERGENCE_ORDER][misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus))[ltsFace];
  auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus))[ltsFace];

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#endif // ACL_DEVICE_OFFLOAD
    for (unsigned i = 0; i < misc::numPaddedPoints; ++i) {
      faultStresses.normalStress[o][i] =
          etaP * (qIMinus[o][6][i] - qIPlus[o][6][i] + qIPlus[o][0][i] * invZp +
                  qIMinus[o][0][i] * invZpNeig);

      faultStresses.traction1[o][i] =
          etaS * (qIMinus[o][7][i] - qIPlus[o][7][i] + qIPlus[o][3][i] * invZs +
                  qIMinus[o][3][i] * invZsNeig);

      faultStresses.traction2[o][i] =
          etaS * (qIMinus[o][8][i] - qIPlus[o][8][i] + qIPlus[o][5][i] * invZs +
                  qIMinus[o][5][i] * invZsNeig);
    }
  }
}
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp end declare target
#endif // ACL_DEVICE_OFFLOAD

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp declare target
#endif // ACL_DEVICE_OFFLOAD
void FrictionSolver::postcomputeImposedStateFromNewStress(const FaultStresses& faultStresses,
                                                          const TractionResults& tractionResults,
                                                          double timeWeights[CONVERGENCE_ORDER],
                                                          unsigned int ltsFace) {
  // this initialization of the kernel could be moved to the initializer
  // set inputParam could be extendent for this (or create own function)
  // the kernel then could be a class attribute and following values are only set once
  //(but be careful of race conditions since this is computed in parallel for each face!!)

  // set imposed state to zero
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#endif // ACL_DEVICE_OFFLOAD
  for (unsigned int i = 0; i < tensor::QInterpolated::size(); i++) {
    imposedStatePlus[ltsFace][i] = static_cast<real>(0.0);
    imposedStateMinus[ltsFace][i] = static_cast<real>(0.0);
  }

  auto invZs = impAndEta[ltsFace].invZs;
  auto invZp = impAndEta[ltsFace].invZp;
  auto invZsNeig = impAndEta[ltsFace].invZsNeig;
  auto invZpNeig = impAndEta[ltsFace].invZpNeig;

  using ImposedStateShapeT = real(*)[misc::numQuantities][misc::numPaddedPoints];
  auto* imposedStateP = (reinterpret_cast<ImposedStateShapeT>(imposedStatePlus))[ltsFace];
  auto* imposedStateM = (reinterpret_cast<ImposedStateShapeT>(imposedStateMinus))[ltsFace];

  using QInterpolatedShapeT =
      real(*)[CONVERGENCE_ORDER][misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus))[ltsFace];
  auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus))[ltsFace];

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    auto weight = timeWeights[o];

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#else
#pragma omp simd
#endif // ACL_DEVICE_OFFLOAD
    for (unsigned i = 0; i < misc::numPaddedPoints; ++i) {
      auto normalStress = faultStresses.normalStress[o][i];
      auto xyTraction = tractionResults.traction1[o][i];
      auto xzTraction = tractionResults.traction2[o][i];

      imposedStateM[0][i] += weight * normalStress;
      imposedStateM[3][i] += weight * xyTraction;
      imposedStateM[5][i] += weight * xzTraction;
      imposedStateM[6][i] +=
          weight * (qIMinus[o][6][i] - invZpNeig * (normalStress - qIMinus[o][0][i]));
      imposedStateM[7][i] +=
          weight * (qIMinus[o][7][i] - invZsNeig * (xyTraction - qIMinus[o][3][i]));
      imposedStateM[8][i] +=
          weight * (qIMinus[o][8][i] - invZsNeig * (xzTraction - qIMinus[o][5][i]));

      imposedStateP[0][i] += weight * normalStress;
      imposedStateP[3][i] += weight * xyTraction;
      imposedStateP[5][i] += weight * xzTraction;
      imposedStateP[6][i] += weight * (qIPlus[o][6][i] + invZp * (normalStress - qIPlus[o][0][i]));
      imposedStateP[7][i] += weight * (qIPlus[o][7][i] + invZs * (xyTraction - qIPlus[o][3][i]));
      imposedStateP[8][i] += weight * (qIPlus[o][8][i] + invZs * (xzTraction - qIPlus[o][5][i]));
    }
  }
}
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp end declare target
#endif // ACL_DEVICE_OFFLOAD

real FrictionSolver::calcSmoothStepIncrement(real currentTime, real dt) {
  real gNuc = calcSmoothStep(currentTime);
  real prevTime = currentTime - dt;
  gNuc = gNuc - calcSmoothStep(prevTime);
  return gNuc;
}

real FrictionSolver::calcSmoothStep(real currentTime) {
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < drParameters.t0) {
    return std::exp(misc::power<2>(currentTime - drParameters.t0) /
                    (currentTime * (currentTime - 2.0 * drParameters.t0)));
  } else {
    return 1.0;
  }
}

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp declare target
#endif // ACL_DEVICE_OFFLOAD
void FrictionSolver::saveRuptureFrontOutput(unsigned int ltsFace) {
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#endif // ACL_DEVICE_OFFLOAD
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    constexpr real ruptureFrontThreshold = 0.001;
    if (ruptureTimePending[ltsFace][pointIndex] &&
        slipRateMagnitude[ltsFace][pointIndex] > ruptureFrontThreshold) {
      ruptureTime[ltsFace][pointIndex] = mFullUpdateTime;
      ruptureTimePending[ltsFace][pointIndex] = false;
    }
  }
}
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp end declare target
#endif // ACL_DEVICE_OFFLOAD

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp declare target
#endif // ACL_DEVICE_OFFLOAD
void FrictionSolver::savePeakSlipRateOutput(unsigned int ltsFace) {

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#endif // ACL_DEVICE_OFFLOAD
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    peakSlipRate[ltsFace][pointIndex] =
        std::max(peakSlipRate[ltsFace][pointIndex], slipRateMagnitude[ltsFace][pointIndex]);
  }
}
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp end declare target
#endif // ACL_DEVICE_OFFLOAD

void FrictionSolver::saveAverageSlipOutput(std::array<real, misc::numPaddedPoints>& tmpSlip,
                                           unsigned int ltsFace) {
  real sumOfTmpSlip = 0;
  if (drParameters.isMagnitudeOutputOn) {
    for (unsigned pointIndex = 0; pointIndex < misc::numberOfBoundaryGaussPoints; pointIndex++) {
      sumOfTmpSlip += tmpSlip[pointIndex];
    }
    averagedSlip[ltsFace] += sumOfTmpSlip / misc::numberOfBoundaryGaussPoints;
  }
}
} // namespace seissol::dr::friction_law
