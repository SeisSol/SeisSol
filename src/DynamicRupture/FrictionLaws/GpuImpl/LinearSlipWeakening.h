#ifndef SEISSOL_GPU_LINEARSLIPWEAKENING_H
#define SEISSOL_GPU_LINEARSLIPWEAKENING_H

#include "DynamicRupture/FrictionLaws/GpuImpl/GpuFrictionSolver.h"

namespace seissol::dr::friction_law::gpu {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <typename Derived>
class LinearSlipWeakeningBase : public GpuFrictionSolver<LinearSlipWeakeningBase<Derived>> {
  public:
  LinearSlipWeakeningBase<Derived>(dr::DRParameters& drParameters)
      : GpuFrictionSolver<LinearSlipWeakeningBase<Derived>>(drParameters){};
  /**
   * critical velocity at which slip rate is considered as being zero for instaneous healing
   */
  static constexpr real u0 = 10e-14;

  real (*dC)[misc::numPaddedPoints];
  real (*muS)[misc::numPaddedPoints];
  real (*muD)[misc::numPaddedPoints];
  real (*cohesion)[misc::numPaddedPoints];

  void updateFrictionAndSlip() {
    auto layerSize{this->currLayerSize};

    #pragma omp target data                       \
    map(to: deltaT[0:CONVERGENCE_ORDER])          \
    device(deviceId)
    {
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
        // computes fault strength, which is the critical value whether active slip exists.
        static_cast<Derived*>(this)->calcStrengthHook(
            this->faultStresses, this->strengthBuffer, timeIndex);
        // computes resulting slip rates, traction and slip dependent on current friction
        // coefficient and strength
        this->calcSlipRateAndTraction(
            this->faultStresses, this->tractionResults, this->strengthBuffer, timeIndex);
        static_cast<Derived*>(this)->calcStateVariableHook(this->stateVariableBuffer, timeIndex);
        this->frictionFunctionHook(this->stateVariableBuffer);
      }
    }
  }

  void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture* dynRup,
                                      real fullUpdateTime) override {
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
  void calcSlipRateAndTraction(FaultStresses* faultStressesPtr,
                               TractionResults* tractionResultsPtr,
                               real (*strengthBuffer)[misc::numPaddedPoints],
                               unsigned int timeIndex) {

    auto layerSize{this->currLayerSize};
    auto* initialStressInFaultCS{this->initialStressInFaultCS};
    auto* impAndEta{this->impAndEta};
    auto* slipRateMagnitude{this->slipRateMagnitude};
    auto* slipRate1{this->slipRate1};
    auto* slipRate2{this->slipRate2};
    auto* traction1{this->traction1};
    auto* traction2{this->traction2};
    auto* slip1{this->slip1};
    auto* slip2{this->slip2};
    auto* deltaT{this->deltaT};

    #pragma omp target teams loop                \
    map(to: deltaT[0:CONVERGENCE_ORDER])         \
    is_device_ptr(faultStressesPtr,              \
                  tractionResultsPtr,            \
                  strengthBuffer,                \
                  initialStressInFaultCS,        \
                  impAndEta,                     \
                  slipRateMagnitude,             \
                  slipRate1,                     \
                  slipRate2,                     \
                  traction1,                     \
                  traction2,                     \
                  slip1,                         \
                  slip2)                         \
    device(deviceId)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      auto& faultStresses = faultStressesPtr[ltsFace];
      auto& tractionResults = tractionResultsPtr[ltsFace];
      auto& strength = strengthBuffer[ltsFace];

      #pragma omp loop bind(parallel)
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
        // calculate absolute value of stress in Y and Z direction
        real totalStress1 = initialStressInFaultCS[ltsFace][pointIndex][3] +
                            faultStresses.traction1[timeIndex][pointIndex];
        real totalStress2 = initialStressInFaultCS[ltsFace][pointIndex][5] +
                            faultStresses.traction2[timeIndex][pointIndex];
        real absoluteShearStress = misc::magnitude(totalStress1, totalStress2);
        // calculate slip rates
        slipRateMagnitude[ltsFace][pointIndex] =
            std::max(static_cast<real>(0.0),
                     (absoluteShearStress - strength[pointIndex]) * impAndEta[ltsFace].invEtaS);
        auto divisor =
            strength[pointIndex] + impAndEta[ltsFace].etaS * slipRateMagnitude[ltsFace][pointIndex];
        slipRate1[ltsFace][pointIndex] =
            slipRateMagnitude[ltsFace][pointIndex] * totalStress1 / divisor;
        slipRate2[ltsFace][pointIndex] =
            slipRateMagnitude[ltsFace][pointIndex] * totalStress2 / divisor;
        // calculate traction
        tractionResults.traction1[timeIndex][pointIndex] =
            faultStresses.traction1[timeIndex][pointIndex] -
            impAndEta[ltsFace].etaS * slipRate1[ltsFace][pointIndex];
        tractionResults.traction2[timeIndex][pointIndex] =
            faultStresses.traction2[timeIndex][pointIndex] -
            impAndEta[ltsFace].etaS * slipRate2[ltsFace][pointIndex];
        traction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
        traction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
        // update directional slip
        slip1[ltsFace][pointIndex] += slipRate1[ltsFace][pointIndex] * deltaT[timeIndex];
        slip2[ltsFace][pointIndex] += slipRate2[ltsFace][pointIndex] * deltaT[timeIndex];
      }
    }
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {
    auto layerSize{this->currLayerSize};
    auto* mu{this->mu};
    auto* muS{this->muS};
    auto* muD{this->muD};

    #pragma omp target teams loop \
    is_device_ptr(stateVariableBuffer, mu, muS, muD) \
    device(deviceId)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      auto& stateVariable = stateVariableBuffer[ltsFace];

      #pragma omp loop bind(parallel)
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
        mu[ltsFace][pointIndex] =
            muS[ltsFace][pointIndex] -
            (muS[ltsFace][pointIndex] - muD[ltsFace][pointIndex]) * stateVariable[pointIndex];
      }
    }
  }

  /**
   * Instantaneous healing option:
   * Reset Mu and Slip, if slipRateMagnitude drops below threshold
   * This function is currently not used, as we miss an appropriate benchmark.
   */
  void instantaneousHealing(unsigned int ltsFace) {
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
  void saveDynamicStressOutput() {
    auto layerSize{this->currLayerSize};
    auto fullUpdateTime{this->mFullUpdateTime};
    auto* dynStressTime{this->dynStressTime};
    auto* dynStressTimePending{this->dynStressTimePending};
    auto* accumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto* dC{this->dC};

    #pragma omp target teams loop           \
    is_device_ptr(dynStressTime,            \
                  dynStressTimePending,     \
                  accumulatedSlipMagnitude, \
                  dC)                       \
    device(deviceId)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      #pragma omp loop bind(parallel)
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
        if (dynStressTimePending[pointIndex] &&
            std::fabs(accumulatedSlipMagnitude[ltsFace][pointIndex]) >= dC[ltsFace][pointIndex]) {
          dynStressTime[ltsFace][pointIndex] = fullUpdateTime;
          dynStressTimePending[ltsFace][pointIndex] = false;
        }
      }
    }
  }

  void preHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {
    static_cast<Derived*>(this)->preHook(stateVariableBuffer);
  }

  void postHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {
    static_cast<Derived*>(this)->postHook(stateVariableBuffer);
  }
};

class LinearSlipWeakeningLaw : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw> {
  public:
  LinearSlipWeakeningLaw(dr::DRParameters& drParameters)
      : LinearSlipWeakeningBase<LinearSlipWeakeningLaw>(drParameters){};

  void preHook(real (*stateVariableBuffer)[misc::numPaddedPoints]){};
  void postHook(real (*stateVariableBuffer)[misc::numPaddedPoints]){};

  void calcStrengthHook(FaultStresses* faultStressesPtr,
                        real (*strengthBuffer)[misc::numPaddedPoints],
                        unsigned int timeIndex);

  void calcStateVariableHook(real (*stateVariable)[misc::numPaddedPoints], unsigned int timeIndex);
};

class LinearSlipWeakeningLawForcedRuptureTime : public LinearSlipWeakeningLaw {
  public:
  LinearSlipWeakeningLawForcedRuptureTime(dr::DRParameters& drParameters)
      : LinearSlipWeakeningLaw(drParameters){};

  real (*forcedRuptureTime)[misc::numPaddedPoints];
  real* tn;

  void preHook(real (*stateVariableBuffer)[misc::numPaddedPoints]);
  void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture* dynRup,
                                      real fullUpdateTime);

  void calcStateVariableHook(real (*stateVariable)[misc::numPaddedPoints], unsigned int timeIndex);
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_LINEARSLIPWEAKENING_H
