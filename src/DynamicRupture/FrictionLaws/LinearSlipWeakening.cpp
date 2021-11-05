#include "LinearSlipWeakening.h"
namespace seissol::dr::friction_law {
void LinearSlipWeakeningLawFL2::calcStrengthHook(std::array<real, numPaddedPoints>& Strength,
                                                 FaultStresses& faultStresses,
                                                 unsigned int timeIndex,
                                                 unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    //-------------------------------------
    // calculate Fault Strength
    // fault strength (Uphoff eq 2.44) with addition cohesion term
    Strength[pointIndex] =
        cohesion[ltsFace][pointIndex] -
        mu[ltsFace][pointIndex] * std::min(initialStressInFaultCS[ltsFace][pointIndex][0] +
                                               faultStresses.NormalStressGP[timeIndex][pointIndex],
                                           static_cast<real>(0.0));
  }
}

void LinearSlipWeakeningLawFL2::calcStateVariableHook(
    std::array<real, numPaddedPoints>& stateVariablePsi,
    std::array<real, numPaddedPoints>& outputSlip,
    dynamicRupture::kernel::resampleParameter& resampleKrnl,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  real resampledSlipRate[numPaddedPoints];
  resampleKrnl.resamplePar = slipRateMagnitude[ltsFace];
  resampleKrnl.resampledPar = resampledSlipRate; // output from execute

  // Resample slip-rate, such that the state (Slip) lies in the same polynomial space as the degrees
  // of freedom resampleMatrix first projects LocSR on the two-dimensional basis on the reference
  // triangle with degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the polynomial
  // at the quadrature points
  resampleKrnl.execute();

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    //-------------------------------------
    // integrate Sliprate To Get Slip = State Variable
    slip[ltsFace][pointIndex] =
        slip[ltsFace][pointIndex] + resampledSlipRate[pointIndex] * deltaT[timeIndex];
    outputSlip[pointIndex] =
        outputSlip[pointIndex] + slipRateMagnitude[ltsFace][pointIndex] * deltaT[timeIndex];

    //-------------------------------------
    // Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
    // actually slip is already the stateVariable for this FL, but to simplify the next equations we
    // divide it here by d_C
    stateVariablePsi[pointIndex] = std::min(
        std::fabs(slip[ltsFace][pointIndex]) / d_c[ltsFace][pointIndex], static_cast<real>(1.0));
  }
}

void LinearSlipWeakeningLawFL16::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                                    seissol::initializers::DynamicRupture* dynRup,
                                                    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  LinearSlipWeakeningLawFL2::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  // maybe change later to const_cast?
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);
  forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime);
  tn = layerData.var(concreteLts->tn);
}

void LinearSlipWeakeningLawFL16::setTimeHook(unsigned int ltsFace) {
  tn[ltsFace] = m_fullUpdateTime;
}

void LinearSlipWeakeningLawFL16::calcStateVariableHook(
    std::array<real, numPaddedPoints>& stateVariablePsi,
    std::array<real, numPaddedPoints>& outputSlip,
    dynamicRupture::kernel::resampleParameter& resampleKrnl,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  LinearSlipWeakeningLawFL2::calcStateVariableHook(
      stateVariablePsi, outputSlip, resampleKrnl, timeIndex, ltsFace);
  tn[ltsFace] += deltaT[timeIndex];

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    real f2 = 0.0;
    if (drParameters.t_0 == 0) {
      if (tn[ltsFace] >= forcedRuptureTime[ltsFace][pointIndex]) {
        f2 = 1.0;
      } else {
        f2 = 0.0;
      }
    } else {
      f2 = std::max(
          static_cast<real>(0.0),
          std::min(static_cast<real>(1.0),
                   (m_fullUpdateTime - forcedRuptureTime[ltsFace][pointIndex]) / drParameters.t_0));
    }
    stateVariablePsi[pointIndex] = std::max(stateVariablePsi[pointIndex], f2);
  }
}

void LinearSlipWeakeningLawBimaterialFL6::calcStrengthHook(
    std::array<real, numPaddedPoints>& strength,
    FaultStresses& faultStresses,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  std::array<real, numPaddedPoints> LocSlipRate;
  std::array<real, numPaddedPoints> sigma;

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    //  modify strength according to prakash clifton
    // literature e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture
    // problems
    LocSlipRate[pointIndex] =
        std::sqrt(slipRateStrike[ltsFace][pointIndex] * slipRateStrike[ltsFace][pointIndex] +
                  slipRateDip[ltsFace][pointIndex] * slipRateDip[ltsFace][pointIndex]);
    sigma[pointIndex] = faultStresses.NormalStressGP[timeIndex][pointIndex] +
                        initialStressInFaultCS[ltsFace][pointIndex][0];
    prak_clif_mod(regularisedStrength[ltsFace][pointIndex],
                  sigma[pointIndex],
                  LocSlipRate[pointIndex],
                  mu[ltsFace][pointIndex],
                  deltaT[timeIndex]);

    // TODO: add this line to make the FL6 actually functional: (this line is also missing in the
    // master branch) strength[pointIndex] = strengthData[ltsFace][pointIndex];
  }
}

void LinearSlipWeakeningLawBimaterialFL6::calcStateVariableHook(
    std::array<real, numPaddedPoints>& stateVariablePsi,
    std::array<real, numPaddedPoints>& outputSlip,
    dynamicRupture::kernel::resampleParameter& resampleKrnl,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    slip[ltsFace][pointIndex] =
        slip[ltsFace][pointIndex] + slipRateMagnitude[ltsFace][pointIndex] * deltaT[timeIndex];
    outputSlip[pointIndex] = slip[ltsFace][pointIndex];

    //-------------------------------------
    // Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
    // actually slip is already the stateVariable for this FL, but to simplify the next equations we
    // divide it here by d_C
    stateVariablePsi[pointIndex] = std::min(
        std::fabs(slip[ltsFace][pointIndex]) / d_c[ltsFace][pointIndex], static_cast<real>(1.0));
  }
}

void LinearSlipWeakeningLawBimaterialFL6::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  LinearSlipWeakeningLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  // maybe change later to const_cast?
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);
  regularisedStrength = layerData.var(concreteLts->regularisedStrength);
}

/*
 * calculates strength
 */
void LinearSlipWeakeningLawBimaterialFL6::prak_clif_mod(
    real& strength, real& sigma, real& LocSlipRate, real& mu, real& dt) {
  real expterm;
  expterm =
      std::exp(-(std::abs(LocSlipRate) + drParameters.v_star) * dt / drParameters.prakash_length);
  strength = strength * expterm - std::max(static_cast<real>(0.0), -mu * sigma) * (expterm - 1.0);
}
} // namespace seissol::dr::friction_law