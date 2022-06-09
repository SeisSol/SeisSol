#include "DynamicRupture/FrictionLaws/GpuImpl/LinearSlipWeakening.h"

namespace seissol::dr::friction_law::gpu {
void LinearSlipWeakeningLaw::calcStrengthHook(FaultStresses* faultStressesPtr,
                                              real (*strengthBuffer)[misc::numPaddedPoints],
                                              unsigned int timeIndex) {

  auto layerSize{this->currLayerSize};
  auto* initialStressInFaultCS{this->initialStressInFaultCS};
  auto* cohesion{this->cohesion};
  auto* mu{this->mu};

  #pragma omp target teams loop              \
  is_device_ptr(faultStressesPtr,            \
                strengthBuffer,              \
                initialStressInFaultCS,      \
                cohesion,                    \
                mu)                          \
  device(deviceId) nowait
  for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
    auto& faultStresses = faultStressesPtr[ltsFace];
    auto& strength = strengthBuffer[ltsFace];

    #pragma omp loop bind(parallel)
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
}

void LinearSlipWeakeningLaw::calcStateVariableHook(
    real (*stateVariableBuffer)[misc::numPaddedPoints], unsigned int timeIndex) {

  auto layerSize{this->currLayerSize};
  auto* accumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
  auto* slipRateMagnitude{this->slipRateMagnitude};
  auto* dC{this->dC};
  auto* deltaT{this->devDeltaT};
  auto* resample{this->resampleMatrix};

  constexpr auto dim0 = misc::dimSize<init::resample, 0>();
  constexpr auto dim1 = misc::dimSize<init::resample, 1>();
  static_assert(dim0 == misc::numPaddedPoints);
  static_assert(dim0 >= dim1);

  #pragma omp target teams loop                 \
  is_device_ptr(stateVariableBuffer,            \
                slipRateMagnitude,              \
                resample,                       \
                accumulatedSlipMagnitude,       \
                dC,                             \
                deltaT)                         \
  device(deviceId) nowait
  for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
    real resampledSlipRate[misc::numPaddedPoints]{};

    // perform matrix vector multiplication
    #pragma omp loop bind(parallel)
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
      real result{0.0};
      for (size_t i{0}; i < dim1; ++i) {
        result += resample[pointIndex + i * dim0] * slipRateMagnitude[ltsFace][i];
      }
      resampledSlipRate[pointIndex] = result;
    }

    auto& stateVariable = stateVariableBuffer[ltsFace];

    #pragma omp loop bind(parallel)
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      // integrate slip rate to get slip = state variable
      accumulatedSlipMagnitude[ltsFace][pointIndex] +=
          resampledSlipRate[pointIndex] * deltaT[timeIndex];

      // Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
      // actually slip is already the stateVariable for this FL, but to simplify the next equations
      // we divide it here by d_C
      stateVariable[pointIndex] = std::min(
          std::fabs(accumulatedSlipMagnitude[ltsFace][pointIndex]) / dC[ltsFace][pointIndex],
          static_cast<real>(1.0));
    }
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
  firstprivate(fullUpdateTime) device(deviceId) nowait
  for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
    tn[ltsFace] = fullUpdateTime;
  }
}

void LinearSlipWeakeningLawForcedRuptureTime::calcStateVariableHook(
    real (*stateVariableBuffer)[misc::numPaddedPoints], unsigned int timeIndex) {
  LinearSlipWeakeningLaw::calcStateVariableHook(stateVariableBuffer, timeIndex);

  auto layerSize{this->currLayerSize};
  auto t0 = drParameters.t0;
  auto* deltaT{this->devDeltaT};
  auto* forcedRuptureTime{this->forcedRuptureTime};
  auto* tn{this->tn};

  #pragma omp target teams loop      \
  is_device_ptr(stateVariableBuffer, \
                forcedRuptureTime,   \
                tn,                  \
                deltaT)              \
  device(deviceId) nowait
  for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
    tn[ltsFace] += deltaT[timeIndex];
    auto& stateVariable = stateVariableBuffer[ltsFace];

    #pragma omp loop bind(parallel)
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      real f2 = 0.0;
      if (t0 == 0) {
        if (tn[ltsFace] >= forcedRuptureTime[ltsFace][pointIndex]) {
          f2 = 1.0;
        } else {
          f2 = 0.0;
        }
      } else {
        f2 = std::max(static_cast<real>(0.0),
                      std::min(static_cast<real>(1.0),
                               // Note: In the fortran implementation on the master branch, this is
                               // m_fullUpdateTime, but this implementation is correct.
                               (tn[ltsFace] - forcedRuptureTime[ltsFace][pointIndex]) / t0));
      }
      stateVariable[pointIndex] = std::max(stateVariable[pointIndex], f2);
    }
  }
}
} // namespace seissol::dr::friction_law::gpu
