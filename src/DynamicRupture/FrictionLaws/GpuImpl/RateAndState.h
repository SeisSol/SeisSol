#ifndef SEISSOL_GPU_RATEANDSTATE_H
#define SEISSOL_GPU_RATEANDSTATE_H

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"
#include <omp.h>

namespace seissol::dr::friction_law::gpu {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived, class TPMethod>
class RateAndStateBase : public BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>> {
  public:
  explicit RateAndStateBase(DRParameters* drParameters)
      : BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>>::BaseFrictionSolver(drParameters),
        tpMethod(TPMethod(drParameters)) {}

  ~RateAndStateBase() {
    if (this->maxClusterSize == 0)
      return;

    omp_free(initialVariables.absoluteShearTraction);
    omp_free(initialVariables.localSlipRate);
    omp_free(initialVariables.normalStress);
    omp_free(initialVariables.stateVarReference);
    omp_free(hasConverged);
  }

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {
    auto* concreteLts = dynamic_cast<seissol::initializers::LTSRateAndState const* const>(dynRup);
    a = layerData.var(concreteLts->rsA);
    sl0 = layerData.var(concreteLts->rsSl0);
    stateVariable = layerData.var(concreteLts->stateVariable);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    tpMethod.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  void allocateAuxiliaryMemory() override {
    FrictionSolverDetails::allocateAuxiliaryMemory();
    if (this->maxClusterSize == 0)
      return;

    {
      using gpPointType = real(*)[misc::numPaddedPoints];
      const size_t requiredNumBytes = misc::numPaddedPoints * this->maxClusterSize * sizeof(real);
      initialVariables.absoluteShearTraction =
          static_cast<gpPointType>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
      initialVariables.localSlipRate =
          static_cast<gpPointType>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
      initialVariables.normalStress =
          static_cast<gpPointType>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
      initialVariables.stateVarReference =
          static_cast<gpPointType>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
    }
    {
      const size_t requiredNumBytes = misc::numPaddedPoints * this->maxClusterSize * sizeof(bool);
      hasConverged = static_cast<bool*>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
    }
  }

  void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture const* const dynRup,
                                      real fullUpdateTime) override {
    auto* concreteLts = dynamic_cast<seissol::initializers::LTSRateAndState const* const>(dynRup);
    this->a = layerData.var(concreteLts->rsA);
    this->sl0 = layerData.var(concreteLts->rsSl0);
    this->stateVariable = layerData.var(concreteLts->stateVariable);
    this->tpMethod.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  #pragma omp declare target
  void updateFrictionAndSlip(unsigned timeIndex) {
    // compute initial slip rate and reference values
    static_cast<Derived*>(this)->calcInitialVariables(timeIndex);

    this->updateStateVariableIterative(timeIndex);
    static_cast<Derived*>(this)->executeIfNotConverged();

    tpMethod.calcFluidPressure(this->initialVariables.normalStress,
                               this->mu,
                               this->initialVariables.localSlipRate,
                               this->deltaT[timeIndex],
                               true);
    updateNormalStress(timeIndex);
    this->calcSlipRateAndTraction(timeIndex);
  }

  void preHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {
    // copy state variable from last time step
    const auto layerSize{this->currLayerSize};

    auto* devLocalStateVariable{this->stateVariable};
    // #pragma omp distribute
    #pragma omp target distribute map(to: devLocalStateVariable[0:layerSize]) map(from: stateVariableBuffer[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        stateVariableBuffer[ltsFace][pointIndex] = devLocalStateVariable[ltsFace][pointIndex];
      }
    }
  }

  void postHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {
    static_cast<Derived*>(this)->resampleStateVar(stateVariableBuffer);
  }

  /**
   * Contains all the variables, which are to be computed initially in each timestep.
   */
  struct InitialVariables {
    real (*absoluteShearTraction)[misc::numPaddedPoints]{nullptr};
    real (*localSlipRate)[misc::numPaddedPoints]{nullptr};
    real (*normalStress)[misc::numPaddedPoints]{nullptr};
    real (*stateVarReference)[misc::numPaddedPoints]{nullptr};
  } initialVariables;

  /*
   * Compute shear stress magnitude, localSlipRate, effective normal stress, reference state
   * variable. Also sets slipRateMagnitude member to reference value.
   */
  void calcInitialVariables(unsigned int timeIndex) {
    const auto layerSize{this->currLayerSize};
    auto* devStateVariableBuffer{this->stateVariableBuffer};
    auto* devFaultStresses{this->faultStresses};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devSlipRate1{this->slipRate1};
    auto* devSlipRate2{this->slipRate2};
    auto* devInitialStressInFaultCS{this->initialStressInFaultCS};

    auto* devAbsoluteShearTraction{this->initialVariables.absoluteShearTraction};
    auto* devLocalSlipRate{this->initialVariables.localSlipRate};
    auto* devStateVarReference{this->initialVariables.stateVarReference};

    updateNormalStress(timeIndex);

    // #pragma omp distribute
    #pragma omp target distribute map(to: devFaultStresses[0:layerSize], devStateVariableBuffer[0:layerSize], devSlipRate1[0:layerSize], devSlipRate2[0:layerSize], devInitialStressInFaultCS[0:layerSize]) map(from: devSlipRateMagnitude[0:layerSize], devAbsoluteShearTraction[0:layerSize], devLocalSlipRate[0:layerSize], devStateVarReference[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        auto& faultStresses = devFaultStresses[ltsFace];

        devStateVarReference[ltsFace][pointIndex] = devStateVariableBuffer[ltsFace][pointIndex];

        const real totalTraction1 = devInitialStressInFaultCS[ltsFace][pointIndex][3] +
                                    faultStresses.traction1[timeIndex][pointIndex];

        const real totalTraction2 = devInitialStressInFaultCS[ltsFace][pointIndex][5] +
                                    faultStresses.traction2[timeIndex][pointIndex];

        devAbsoluteShearTraction[ltsFace][pointIndex] =
            misc::magnitude(totalTraction1, totalTraction2);
        auto localSlipRateMagnitude =
            misc::magnitude(devSlipRate1[ltsFace][pointIndex], devSlipRate2[ltsFace][pointIndex]);

        localSlipRateMagnitude = std::max(rs::almostZero(), localSlipRateMagnitude);
        devSlipRateMagnitude[ltsFace][pointIndex] = localSlipRateMagnitude;
        devLocalSlipRate[ltsFace][pointIndex] = localSlipRateMagnitude;
      }
    }
  }

  void updateStateVariableIterative(unsigned timeIndex) {
    const auto layerSize{this->currLayerSize};
    auto* devHasConverged{this->hasConverged};
    auto* devLocalSlipRate{this->initialVariables.localSlipRate};
    auto* devStateVariableBuffer{this->stateVariableBuffer};
    auto* devNormalStress{this->initialVariables.normalStress};
    auto* devAbsoluteShearStress{this->initialVariables.absoluteShearTraction};
    auto* devMu{this->mu};

    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devImpAndEta{this->impAndEta};
    auto solverSettings{this->settings};

    auto details = static_cast<Derived*>(this)->getCurrentLtsLayerDetails();
    for (unsigned j = 0; j < this->settings.numberStateVariableUpdates; j++) {

      const auto dt{this->deltaT[timeIndex]};
      static_cast<Derived*>(this)->updateStateVariable(dt);
      this->tpMethod.calcFluidPressure(devNormalStress, devMu, devLocalSlipRate, dt, false);
      updateNormalStress(timeIndex);

      // #pragma omp distribute
      #pragma omp target distribute map(to: devStateVariableBuffer[0:layerSize], devNormalStress[0:layerSize], devAbsoluteShearStress[0:layerSize], devImpAndEta[0:layerSize]) map(tofrom: devSlipRateMagnitude[0:layerSize]) map(from: devHasConverged[0:layerSize], devLocalSlipRate[0:layerSize], devMu[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        bool hasConvergedAllPoints = true;
        #pragma omp parallel for schedule(static, 1) reduction(&&:hasConvergedAllPoints)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        const auto localStateVariable = devStateVariableBuffer[ltsFace][pointIndex];
        const auto normalStress = devNormalStress[ltsFace][pointIndex];
        const auto absoluteShearStress = devAbsoluteShearStress[ltsFace][pointIndex];
        const auto localSlipRateMagnitude = devSlipRateMagnitude[ltsFace][pointIndex];
        const auto localImpAndEta = devImpAndEta[ltsFace];

        real slipRateTest{};
        bool hasConvergedLocal = RateAndStateBase::invertSlipRateIterative(slipRateTest,
                                                                            localStateVariable,
                                                                            normalStress,
                                                                            absoluteShearStress,
                                                                            localSlipRateMagnitude,
                                                                            localImpAndEta.invEtaS,
                                                                            details,
                                                                            solverSettings,
                                                                            ltsFace,
                                                                            pointIndex);

        hasConvergedAllPoints &= hasConvergedLocal;

          devLocalSlipRate[ltsFace][pointIndex] =
              0.5 * (localSlipRateMagnitude + std::fabs(slipRateTest));
          devSlipRateMagnitude[ltsFace][pointIndex] = std::fabs(slipRateTest);

          devMu[ltsFace][pointIndex] = Derived::updateMu(
              localSlipRateMagnitude, localStateVariable, details, ltsFace, pointIndex);
        }
        devHasConverged[ltsFace] = hasConvergedAllPoints;
      }
    }
  }

  void calcSlipRateAndTraction(unsigned timeIndex) {
    const auto layerSize{this->currLayerSize};
    auto* devStateVarReference{this->initialVariables.stateVarReference};
    auto* devLocalSlipRate{this->initialVariables.localSlipRate};
    auto* devStateVariableBuffer{this->stateVariableBuffer}; // localStateVariable
    auto* devNormalStress{this->initialVariables.normalStress};
    auto* devAbsoluteTraction{this->initialVariables.absoluteShearTraction};
    auto* devFaultStresses{this->faultStresses};
    auto* devTractionResults{this->tractionResults};

    auto* devMu{this->mu};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
    auto* devTraction1{this->traction1};
    auto* devTraction2{this->traction2};
    auto* devSlipRate1{this->slipRate1};
    auto* devSlipRate2{this->slipRate2};
    auto* devSlip1{this->slip1};
    auto* devSlip2{this->slip2};
    auto* devImpAndEta{this->impAndEta};
    auto* devAccumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto deltaTime{this->deltaT[timeIndex]};

    auto details = static_cast<Derived*>(this)->getCurrentLtsLayerDetails();

    static_cast<Derived*>(this)->updateStateVariable(this->deltaT[timeIndex]);

    // #pragma omp distribute
    #pragma omp target distribute map(to: devStateVariableBuffer[0:layerSize], devSlipRateMagnitude[0:layerSize], devNormalStress[0:layerSize], devAbsoluteTraction[0:layerSize], devFaultStresses[0:layerSize], devInitialStressInFaultCS[0:layerSize], devImpAndEta[0:layerSize]) map(tofrom: devMu[0:layerSize], devAccumulatedSlipMagnitude[0:layerSize]) map(from: devTraction1[0:layerSize], devTraction2[0:layerSize], devSlip1[0:layerSize], devSlip2[0:layerSize], devTractionResults[0:layerSize], devSlipRate1[0:layerSize], devSlipRate2[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {

        const auto localStateVariable = devStateVariableBuffer[ltsFace][pointIndex];
        const auto slipRateMagnitude = devSlipRateMagnitude[ltsFace][pointIndex];
        devMu[ltsFace][pointIndex] =
            Derived::updateMu(slipRateMagnitude, localStateVariable, details, ltsFace, pointIndex);
        const real strength = -devMu[ltsFace][pointIndex] * devNormalStress[ltsFace][pointIndex];

        const auto initialStressInFaultCS = devInitialStressInFaultCS[ltsFace][pointIndex];
        const auto savedTraction1 = devFaultStresses[ltsFace].traction1[timeIndex][pointIndex];
        const auto savedTraction2 = devFaultStresses[ltsFace].traction2[timeIndex][pointIndex];

        // calculate absolute value of stress in Y and Z direction
        devFaultStresses[ltsFace].traction2[timeIndex][pointIndex];
        const real totalTraction1 = initialStressInFaultCS[3] + savedTraction1;
        const real totalTraction2 = initialStressInFaultCS[5] + savedTraction2;

        // update stress change
        const auto traction1 =
            (totalTraction1 / devAbsoluteTraction[ltsFace][pointIndex]) * strength -
            initialStressInFaultCS[3];
        const auto traction2 =
            (totalTraction2 / devAbsoluteTraction[ltsFace][pointIndex]) * strength -
            initialStressInFaultCS[5];

        // Compute slip
        devAccumulatedSlipMagnitude[ltsFace][pointIndex] += slipRateMagnitude * deltaTime;

        // Update slip rate
        const auto invEtaS = devImpAndEta[ltsFace].invEtaS;
        auto slipRate1 = -invEtaS * (traction1 - savedTraction1);
        auto slipRate2 = -invEtaS * (traction2 - savedTraction2);

        const real locSlipRateMagnitude = misc::magnitude(slipRate1, slipRate2);

        if (locSlipRateMagnitude != 0.0) {
          slipRate1 *= slipRateMagnitude / locSlipRateMagnitude;
          slipRate2 *= slipRateMagnitude / locSlipRateMagnitude;
        }

        // Save traction for flux computation
        devTraction1[ltsFace][pointIndex] = traction1;
        devTraction2[ltsFace][pointIndex] = traction2;

        // update directional slip
        devSlip1[ltsFace][pointIndex] += slipRate1 * deltaTime;
        devSlip2[ltsFace][pointIndex] += slipRate2 * deltaTime;

        // update traction
        devTractionResults[ltsFace].traction1[timeIndex][pointIndex] = traction1;
        devTractionResults[ltsFace].traction2[timeIndex][pointIndex] = traction2;

        // update slip rate
        devSlipRate1[ltsFace][pointIndex] = slipRate1;
        devSlipRate2[ltsFace][pointIndex] = slipRate2;
      }
    }
  }

  void saveDynamicStressOutput() {
    const auto layerSize{this->currLayerSize};
    auto fullUpdateTime{this->mFullUpdateTime};
    auto muW{this->drParameters->muW};
    auto rsF0{this->drParameters->rsF0};

    auto* devDynStressTime{this->dynStressTime};
    auto* devDynStressTimePending{this->dynStressTimePending};
    auto* devRuptureTime{this->ruptureTime};
    auto* devMu{this->mu};

    // #pragma omp distribute
    #pragma omp target distribute map(to: devMu[0:layerSize], devRuptureTime[0:layerSize]) map(tofrom: devDynStressTimePending[0:layerSize]) map(from: devDynStressTime[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {

        const auto localRuptureTime = devRuptureTime[ltsFace][pointIndex];
        if (localRuptureTime > 0.0 && localRuptureTime <= fullUpdateTime &&
            devDynStressTimePending[ltsFace][pointIndex] &&
            devMu[ltsFace][pointIndex] <= (muW + 0.05 * (rsF0 - muW))) {
          devDynStressTime[ltsFace][pointIndex] = fullUpdateTime;
          devDynStressTimePending[ltsFace][pointIndex] = false;
        }
      }
    }
  }

  template <typename DetailsType>
  static bool invertSlipRateIterative(real& slipRateTest,
                                      real localStateVariable,
                                      real normalStress,
                                      real absoluteShearStress,
                                      real slipRateMagnitude,
                                      real invEtaS,
                                      DetailsType details,
                                      rs::Settings solverSettings,
                                      int ltsFace,
                                      int pointIndex) {

    // Note that we need double precision here, since single precision led to NaNs.
    double muF{0.0}, dMuF{0.0}, g{0.0}, dG{0.0};
    slipRateTest = slipRateMagnitude;

    for (unsigned i = 0; i < solverSettings.maxNumberSlipRateUpdates; i++) {
      bool converged;
        muF = Derived::updateMu(slipRateTest, localStateVariable, details, ltsFace, pointIndex);
        dMuF = Derived::updateMuDerivative(
            slipRateTest, localStateVariable, details, ltsFace, pointIndex);

        g = -invEtaS * (std::fabs(normalStress) * muF - absoluteShearStress) - slipRateTest;
        converged = std::fabs(g) < solverSettings.newtonTolerance;

      if (converged) { return true; }

        dG = -invEtaS * (std::fabs(normalStress) * dMuF) - 1.0;
        slipRateTest =
            std::max(friction_law::rs::almostZero(), static_cast<real>(slipRateTest - (g / dG)));
    }
    return false;
  }

  void updateNormalStress(size_t timeIndex) {
    const auto layerSize{this->currLayerSize};
    auto* devFaultStresses{this->faultStresses};
    auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
    auto* devNormalStress{this->initialVariables.normalStress};

    auto tpCurrentLayerDetails = tpMethod.getCurrentLayerDetails();

    // #pragma omp distribute
    #pragma omp target distribute map(to: devFaultStresses[0:layerSize], devInitialStressInFaultCS[0:layerSize], tpCurrentLayerDetails) map(from: devNormalStress[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        auto& faultStresses = devFaultStresses[ltsFace];

        devNormalStress[ltsFace][pointIndex] =
            std::min(static_cast<real>(0.0),
                     faultStresses.normalStress[timeIndex][pointIndex] +
                         devInitialStressInFaultCS[ltsFace][pointIndex][0] -
                         TPMethod::getFluidPressure(tpCurrentLayerDetails, ltsFace, pointIndex));
      }
    }
  }

  protected:
  real (*a)[misc::numPaddedPoints];
  real (*sl0)[misc::numPaddedPoints];
  real (*stateVariable)[misc::numPaddedPoints];
  bool* hasConverged{};

  TPMethod tpMethod;
  rs::Settings settings{};
  #pragma omp end declare target
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_RATEANDSTATE_H
