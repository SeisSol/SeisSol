#include "LinearSlipWeakening.h"
namespace seissol::dr::friction_law {
void LinearSlipWeakeningLaw::calcStrengthHook(FaultStresses& faultStresses,
                                              std::array<real, numPaddedPoints>& strength,
                                              unsigned int timeIndex,
                                              unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    //-------------------------------------
    // calculate Fault Strength
    // fault strength (Uphoff eq 2.44) with addition cohesion term
    real totalNormalStress = initialStressInFaultCS[ltsFace][pointIndex][0] +
                             faultStresses.NormalStressGP[timeIndex][pointIndex];
    strength[pointIndex] =
        -cohesion[ltsFace][pointIndex] -
        mu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<real>(0.0));
  }
}

void LinearSlipWeakeningLaw::calcStateVariableHook(std::array<real, numPaddedPoints>& stateVariable,
                                                   unsigned int timeIndex,
                                                   unsigned int ltsFace) {
  real resampledSlipRate[numPaddedPoints];
  dynamicRupture::kernel::resampleParameter resampleKrnl;
  resampleKrnl.resampleM = init::resample::Values;
  resampleKrnl.resamplePar = slipRateMagnitude[ltsFace];
  resampleKrnl.resampledPar = resampledSlipRate;

  // Resample slip-rate, such that the state increment (slip) lies in the same polynomial space as
  // the degrees of freedom resampleMatrix first projects LocSR on the two-dimensional basis on the
  // reference triangle with degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the
  // polynomial at the quadrature points
  resampleKrnl.execute();

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    // integrate slip rate to get slip = state variable
    accumulatedSlipMagnitude[ltsFace][pointIndex] =
        accumulatedSlipMagnitude[ltsFace][pointIndex] +
        resampledSlipRate[pointIndex] * deltaT[timeIndex];

    // Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
    // actually slip is already the stateVariable for this FL, but to simplify the next equations we
    // divide it here by d_C
    stateVariable[pointIndex] = std::min(std::fabs(accumulatedSlipMagnitude[ltsFace][pointIndex]) /
                                             d_c[ltsFace][pointIndex],
                                         static_cast<real>(1.0));
  }
}

void LinearSlipWeakeningLaw::preHook(std::array<real, numPaddedPoints>& stateVariableBuffer,
                                     unsigned int ltsFace) {}
void LinearSlipWeakeningLaw::postHook(std::array<real, numPaddedPoints>& stateVariableBuffer,
                                      unsigned int ltsFace) {}

void LinearSlipWeakeningLawForcedRuptureTime::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  LinearSlipWeakeningLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  // maybe change later to const_cast?
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);
  forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime);
  tn = layerData.var(concreteLts->tn);
}

void LinearSlipWeakeningLawForcedRuptureTime::preHook(
    std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace) {
  tn[ltsFace] = m_fullUpdateTime;
}

void LinearSlipWeakeningLawForcedRuptureTime::calcStateVariableHook(
    std::array<real, numPaddedPoints>& stateVariable,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  LinearSlipWeakeningLaw::calcStateVariableHook(stateVariable, timeIndex, ltsFace);
  tn[ltsFace] += deltaT[timeIndex];

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    real f2 = 0.0;
    if (drParameters.t0 == 0) {
      if (tn[ltsFace] >= forcedRuptureTime[ltsFace][pointIndex]) {
        f2 = 1.0;
      } else {
        f2 = 0.0;
      }
    } else {
      f2 = std::max(
          static_cast<real>(0.0),
          std::min(static_cast<real>(1.0),
                   // Note: In the fortran implementation on the master branch, this is
                   // m_fullUpdateTime, but this implementation is correct.
                   (tn[ltsFace] - forcedRuptureTime[ltsFace][pointIndex]) / drParameters.t0));
    }
    stateVariable[pointIndex] = std::max(stateVariable[pointIndex], f2);
  }
}

void LinearSlipWeakeningLawBimaterial::calcStrengthHook(FaultStresses& faultStresses,
                                                        std::array<real, numPaddedPoints>& strength,
                                                        unsigned int timeIndex,
                                                        unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    //  modify strength according to prakash clifton
    // literature e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture
    // problems
    real localSlipRate =
        misc::magnitude(slipRate1[ltsFace][pointIndex], slipRate2[ltsFace][pointIndex]);
    real sigma = faultStresses.NormalStressGP[timeIndex][pointIndex] +
                 initialStressInFaultCS[ltsFace][pointIndex][0];
    prak_clif_mod(regularisedStrength[ltsFace][pointIndex],
                  sigma,
                  localSlipRate,
                  mu[ltsFace][pointIndex],
                  deltaT[timeIndex]);

    // TODO: add this line to make the FL6 actually functional: (this line is also missing in the
    // master branch)
    strength[pointIndex] = regularisedStrength[ltsFace][pointIndex];
  }
}

void LinearSlipWeakeningLawBimaterial::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  LinearSlipWeakeningBase::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  // maybe change later to const_cast?
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);
  regularisedStrength = layerData.var(concreteLts->regularisedStrength);
}

/*
 * calculates strength
 */
void LinearSlipWeakeningLawBimaterial::prak_clif_mod(
    real& strength, real& sigma, real& localSlipRate, real& mu, real& dt) {
  real expterm =
      std::exp(-(std::abs(localSlipRate) + drParameters.vStar) * dt / drParameters.prakashLength);
  strength = strength * expterm - std::max(static_cast<real>(0.0), -mu * sigma) * (expterm - 1.0);
}
} // namespace seissol::dr::friction_law