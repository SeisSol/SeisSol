#include "DynamicRupture/FrictionLaws/GpuImpl/LinearSlipWeakening.h"

namespace seissol::dr::friction_law::gpu {
void LinearSlipWeakeningBase::updateFrictionAndSlip(
    FaultStresses& faultStresses,
    TractionResults& tractionResults,
    std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
    std::array<real, misc::numPaddedPoints>& strengthBuffer,
    unsigned int ltsFace,
    unsigned int timeIndex) {
  // computes fault strength, which is the critical value whether active slip exists.
  this->calcStrengthHook(faultStresses, strengthBuffer, timeIndex, ltsFace);
  // computes resulting slip rates, traction and slip dependent on current friction
  // coefficient and strength
  this->calcSlipRateAndTraction(faultStresses, tractionResults, strengthBuffer, timeIndex, ltsFace);
  this->calcStateVariableHook(stateVariableBuffer, timeIndex, ltsFace);
  this->frictionFunctionHook(stateVariableBuffer, ltsFace);
}

void LinearSlipWeakeningBase::copySpecificLtsDataTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
  this->dC = layerData.var(concreteLts->dC);
  this->muS = layerData.var(concreteLts->muS);
  this->muD = layerData.var(concreteLts->muD);
  this->cohesion = layerData.var(concreteLts->cohesion);
}

/**
 *  compute the slip rate and the traction from the fault strength and fault stresses
 *  also updates the directional slip1 and slip2
 */
void LinearSlipWeakeningBase::calcSlipRateAndTraction(
    FaultStresses& faultStresses,
    TractionResults& tractionResults,
    std::array<real, misc::numPaddedPoints>& strength,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    // calculate absolute value of stress in Y and Z direction
    real totalStressXY = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                         faultStresses.xyStress[timeIndex][pointIndex];
    real totalStressXZ = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                         faultStresses.xzStress[timeIndex][pointIndex];
    real absoluteShearStress = misc::magnitude(totalStressXY, totalStressXZ);
    // calculate slip rates
    this->slipRateMagnitude[ltsFace][pointIndex] =
        std::max(static_cast<real>(0.0),
                 (absoluteShearStress - strength[pointIndex]) * this->impAndEta[ltsFace].invEtaS);
    auto divisor = strength[pointIndex] +
                   this->impAndEta[ltsFace].etaS * this->slipRateMagnitude[ltsFace][pointIndex];
    this->slipRate1[ltsFace][pointIndex] =
        this->slipRateMagnitude[ltsFace][pointIndex] * totalStressXY / divisor;
    this->slipRate2[ltsFace][pointIndex] =
        this->slipRateMagnitude[ltsFace][pointIndex] * totalStressXZ / divisor;
    // calculate traction
    tractionResults.xyTraction[timeIndex][pointIndex] =
        faultStresses.xyStress[timeIndex][pointIndex] -
        this->impAndEta[ltsFace].etaS * this->slipRate1[ltsFace][pointIndex];
    tractionResults.xzTraction[timeIndex][pointIndex] =
        faultStresses.xzStress[timeIndex][pointIndex] -
        this->impAndEta[ltsFace].etaS * this->slipRate2[ltsFace][pointIndex];
    this->tractionXY[ltsFace][pointIndex] = tractionResults.xyTraction[timeIndex][pointIndex];
    this->tractionXZ[ltsFace][pointIndex] = tractionResults.xzTraction[timeIndex][pointIndex];
    // update directional slip
    this->slip1[ltsFace][pointIndex] +=
        this->slipRate1[ltsFace][pointIndex] * this->deltaT[timeIndex];
    this->slip2[ltsFace][pointIndex] +=
        this->slipRate2[ltsFace][pointIndex] * this->deltaT[timeIndex];
  }
}
/**
 * evaluate friction law: updated mu -> friction law
 * for example see Carsten Uphoff's thesis: Eq. 2.45
 */
void LinearSlipWeakeningBase::frictionFunctionHook(
    std::array<real, misc::numPaddedPoints>& stateVariable, unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    this->mu[ltsFace][pointIndex] =
        muS[ltsFace][pointIndex] -
        (muS[ltsFace][pointIndex] - muD[ltsFace][pointIndex]) * stateVariable[pointIndex];
  }
}

/**
 * Instantaneous healing option:
 * Reset Mu and Slip, if slipRateMagnitude drops below threshold
 * This function is currently not used, as we miss an appropriate benchmark.
 */
void LinearSlipWeakeningBase::instantaneousHealing(unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    if (this->slipRateMagnitude[ltsFace][pointIndex] < u0) {
      this->mu[ltsFace][pointIndex] = muS[ltsFace][pointIndex];
      this->accumulatedSlipMagnitude[ltsFace][pointIndex] = 0.0;
    }
  }
}

/*
 * output time when shear stress is equal to the dynamic stress after rupture arrived
 * currently only for linear slip weakening
 */
void LinearSlipWeakeningBase::saveDynamicStressOutput(unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    if (this->dynStressTimePending[pointIndex] &&
        std::fabs(this->accumulatedSlipMagnitude[ltsFace][pointIndex]) >= dC[ltsFace][pointIndex]) {
      this->dynStressTime[ltsFace][pointIndex] = this->mFullUpdateTime;
      this->dynStressTimePending[ltsFace][pointIndex] = false;
    }
  }
}

//--------------------------------------------------------------------------

void LinearSlipWeakeningLaw::calcStrengthHook(FaultStresses& faultStresses,
                                              std::array<real, misc::numPaddedPoints>& strength,
                                              unsigned int timeIndex,
                                              unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    //-------------------------------------
    // calculate Fault Strength
    // fault strength (Uphoff eq 2.44) with addition cohesion term
    real totalNormalStress = initialStressInFaultCS[ltsFace][pointIndex][0] +
                             faultStresses.normalStress[timeIndex][pointIndex];
    strength[pointIndex] =
        -cohesion[ltsFace][pointIndex] -
        mu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<real>(0.0));
  }
}

void LinearSlipWeakeningLaw::calcStateVariableHook(
    std::array<real, misc::numPaddedPoints>& stateVariable,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  alignas(ALIGNMENT) real resampledSlipRate[misc::numPaddedPoints]{};
  dynamicRupture::kernel::resampleParameter resampleKrnl;
  resampleKrnl.resample = init::resample::Values;
  resampleKrnl.originalQ = slipRateMagnitude[ltsFace];
  resampleKrnl.resampledQ = resampledSlipRate;

  // Resample slip-rate, such that the state increment (slip) lies in the same polynomial space as
  // the degrees of freedom resampleMatrix first projects LocSR on the two-dimensional basis on the
  // reference triangle with degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the
  // polynomial at the quadrature points
  resampleKrnl.execute();

  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    // integrate slip rate to get slip = state variable
    accumulatedSlipMagnitude[ltsFace][pointIndex] +=
        resampledSlipRate[pointIndex] * deltaT[timeIndex];

    // Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
    // actually slip is already the stateVariable for this FL, but to simplify the next equations we
    // divide it here by d_C
    stateVariable[pointIndex] =
        std::min(std::fabs(accumulatedSlipMagnitude[ltsFace][pointIndex]) / dC[ltsFace][pointIndex],
                 static_cast<real>(1.0));
  }
}

void LinearSlipWeakeningLawForcedRuptureTime::copySpecificLtsDataTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  LinearSlipWeakeningLaw::copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);
  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);
  forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime);
  tn = layerData.var(concreteLts->tn);
}

void LinearSlipWeakeningLawForcedRuptureTime::preHook(
    std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace) {
  tn[ltsFace] = mFullUpdateTime;
}

void LinearSlipWeakeningLawForcedRuptureTime::calcStateVariableHook(
    std::array<real, misc::numPaddedPoints>& stateVariable,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  LinearSlipWeakeningLaw::calcStateVariableHook(stateVariable, timeIndex, ltsFace);
  tn[ltsFace] += deltaT[timeIndex];

  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
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
} // namespace seissol::dr::friction_law::gpu