#include "DynamicRupture/FrictionLaws/GpuImpl/LinearSlipWeakening.h"

namespace seissol::dr::friction_law::gpu {
void LinearSlipWeakeningLaw::calcStrengthHook(FaultStresses& faultStresses,
                                              real (*strength)[misc::numPaddedPoints],
                                              unsigned int timeIndex,
                                              unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    //-------------------------------------
    // calculate Fault Strength
    // fault strength (Uphoff eq 2.44) with addition cohesion term
    real totalNormalStress = initialStressInFaultCS[ltsFace][pointIndex][0] +
                             faultStresses.normalStress[timeIndex][pointIndex];
    strength[ltsFace][pointIndex] =
        -cohesion[ltsFace][pointIndex] -
        mu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<real>(0.0));
  }
}

void LinearSlipWeakeningLaw::calcStateVariableHook(real (*stateVariable)[misc::numPaddedPoints],
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
    stateVariable[ltsFace][pointIndex] =
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
    real (*stateVariableBuffer)[misc::numPaddedPoints]) {
  auto layerSize{this->currLayerSize};
  auto fullUpdateTime{this->mFullUpdateTime};
  auto* tn{this->tn};

  #pragma omp target teams loop is_device_ptr(tn) \
  firstprivate(fullUpdateTime) device(diviceId)
  for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
    tn[ltsFace] = mFullUpdateTime;
  }
}

void LinearSlipWeakeningLawForcedRuptureTime::calcStateVariableHook(
    real (*stateVariable)[misc::numPaddedPoints], unsigned int timeIndex, unsigned int ltsFace) {
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
    stateVariable[ltsFace][pointIndex] = std::max(stateVariable[ltsFace][pointIndex], f2);
  }
}
} // namespace seissol::dr::friction_law::gpu
