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

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             TractionResults& tractionResults,
                             real (*stateVariableBuffer)[misc::numPaddedPoints],
                             real (*strengthBuffer)[misc::numPaddedPoints],
                             unsigned int ltsFace,
                             unsigned int timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    static_cast<Derived*>(this)->calcStrengthHook(
        faultStresses, strengthBuffer, timeIndex, ltsFace);
    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    this->calcSlipRateAndTraction(
        faultStresses, tractionResults, strengthBuffer, timeIndex, ltsFace);
    static_cast<Derived*>(this)->calcStateVariableHook(stateVariableBuffer, timeIndex, ltsFace);
    this->frictionFunctionHook(stateVariableBuffer, ltsFace);
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
  void calcSlipRateAndTraction(FaultStresses& faultStresses,
                               TractionResults& tractionResults,
                               real (*strength)[misc::numPaddedPoints],
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      real totalStressXY = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                           faultStresses.traction1[timeIndex][pointIndex];
      real totalStressXZ = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                           faultStresses.traction2[timeIndex][pointIndex];
      real absoluteShearStress = misc::magnitude(totalStressXY, totalStressXZ);
      // calculate slip rates
      this->slipRateMagnitude[ltsFace][pointIndex] = std::max(
          static_cast<real>(0.0),
          (absoluteShearStress - strength[ltsFace][pointIndex]) * this->impAndEta[ltsFace].invEtaS);
      auto divisor = strength[ltsFace][pointIndex] +
                     this->impAndEta[ltsFace].etaS * this->slipRateMagnitude[ltsFace][pointIndex];
      this->slipRate1[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalStressXY / divisor;
      this->slipRate2[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalStressXZ / divisor;
      // calculate traction
      tractionResults.traction1[timeIndex][pointIndex] =
          faultStresses.traction1[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * this->slipRate1[ltsFace][pointIndex];
      tractionResults.traction2[timeIndex][pointIndex] =
          faultStresses.traction2[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * this->slipRate2[ltsFace][pointIndex];
      this->traction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
      this->traction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
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
  void frictionFunctionHook(real (*stateVariable)[misc::numPaddedPoints], unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      this->mu[ltsFace][pointIndex] =
          muS[ltsFace][pointIndex] - (muS[ltsFace][pointIndex] - muD[ltsFace][pointIndex]) *
                                         stateVariable[ltsFace][pointIndex];
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
    auto* dC = this->dC;
    auto* dynStressTime = this->dynStressTime;
    auto* dynStressTimePending = this->dynStressTimePending;
    auto* accumulatedSlipMagnitude = this->accumulatedSlipMagnitude;

    #pragma omp target teams loop           \
    is_device_ptr(dynStressTimePending,     \
                  accumulatedSlipMagnitude, \
                  dynStressTime,            \
                  dC)                       \
    device(diviceId)
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

  void calcStrengthHook(FaultStresses& faultStresses,
                        real (*strength)[misc::numPaddedPoints],
                        unsigned int timeIndex,
                        unsigned int ltsFace);

  void calcStateVariableHook(real (*stateVariable)[misc::numPaddedPoints],
                             unsigned int timeIndex,
                             unsigned int ltsFace);
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

  void calcStateVariableHook(real (*stateVariable)[misc::numPaddedPoints],
                             unsigned int timeIndex,
                             unsigned int ltsFace);
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_LINEARSLIPWEAKENING_H
