#ifndef SEISSOL_GPU_RATEANDSTATE_H
#define SEISSOL_GPU_RATEANDSTATE_H

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"

namespace seissol::dr::friction_law::gpu {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived, class TPMethod>
class RateAndStateBase : public BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>> {
  public:
  explicit RateAndStateBase(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>>::BaseFrictionSolver(drParameters),
        tpMethod(TPMethod(drParameters)) {}

  ~RateAndStateBase() {
    if (this->maxClusterSize == 0)
      return;

    sycl::free(initialVariables.absoluteShearTraction, this->queue);
    sycl::free(initialVariables.localSlipRate, this->queue);
    sycl::free(initialVariables.normalStress, this->queue);
    sycl::free(initialVariables.stateVarReference, this->queue);
    sycl::free(hasConverged, this->queue);
    this->queue.wait_and_throw();
  }

  void allocateAuxiliaryMemory() override {
    FrictionSolverDetails::allocateAuxiliaryMemory();
    if (this->maxClusterSize == 0)
      return;

    {
      using gpPointType = real(*)[misc::NumPaddedPoints];
      const size_t requiredNumBytes = misc::NumPaddedPoints * this->maxClusterSize * sizeof(real);
      initialVariables.absoluteShearTraction =
          static_cast<gpPointType>(sycl::malloc_device(requiredNumBytes, this->queue));
      initialVariables.localSlipRate =
          static_cast<gpPointType>(sycl::malloc_device(requiredNumBytes, this->queue));
      initialVariables.normalStress =
          static_cast<gpPointType>(sycl::malloc_device(requiredNumBytes, this->queue));
      initialVariables.stateVarReference =
          static_cast<gpPointType>(sycl::malloc_device(requiredNumBytes, this->queue));
    }
    {
      const size_t requiredNumBytes = misc::NumPaddedPoints * this->maxClusterSize * sizeof(bool);
      hasConverged = static_cast<bool*>(sycl::malloc_device(requiredNumBytes, this->queue));
    }
  }

  void copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                      const seissol::initializer::DynamicRupture* const dynRup,
                                      real fullUpdateTime) override {
    auto* concreteLts = dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup);
    this->a = layerData.var(concreteLts->rsA, seissol::initializer::AllocationPlace::Device);
    this->sl0 = layerData.var(concreteLts->rsSl0, seissol::initializer::AllocationPlace::Device);
    this->stateVariable =
        layerData.var(concreteLts->stateVariable, seissol::initializer::AllocationPlace::Device);
    this->tpMethod.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

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

  void preHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {
    // copy state variable from last time step

    auto* devLocalStateVariable{this->stateVariable};
    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);
        stateVariableBuffer[ltsFace][pointIndex] = devLocalStateVariable[ltsFace][pointIndex];
      });
    });
  }

  void postHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {
    static_cast<Derived*>(this)->resampleStateVar(stateVariableBuffer);
  }

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {
    auto* concreteLts = dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup);
    a = layerData.var(concreteLts->rsA, seissol::initializer::AllocationPlace::Device);
    sl0 = layerData.var(concreteLts->rsSl0, seissol::initializer::AllocationPlace::Device);
    stateVariable =
        layerData.var(concreteLts->stateVariable, seissol::initializer::AllocationPlace::Device);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    tpMethod.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  /**
   * Contains all the variables, which are to be computed initially in each timestep.
   */
  struct InitialVariables {
    real (*absoluteShearTraction)[misc::NumPaddedPoints]{nullptr};
    real (*localSlipRate)[misc::NumPaddedPoints]{nullptr};
    real (*normalStress)[misc::NumPaddedPoints]{nullptr};
    real (*stateVarReference)[misc::NumPaddedPoints]{nullptr};
  } initialVariables;

  /*
   * Compute shear stress magnitude, localSlipRate, effective normal stress, reference state
   * variable. Also sets slipRateMagnitude member to reference value.
   */
  void calcInitialVariables(unsigned int timeIndex) {
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

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);
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
      });
    });
  }

  void updateStateVariableIterative(unsigned timeIndex) {
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
    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    for (unsigned j = 0; j < this->settings.numberStateVariableUpdates; j++) {

      const auto dt{this->deltaT[timeIndex]};
      static_cast<Derived*>(this)->updateStateVariable(dt);
      this->tpMethod.calcFluidPressure(devNormalStress, devMu, devLocalSlipRate, dt, false);
      updateNormalStress(timeIndex);

      this->queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          const auto ltsFace = item.get_group().get_group_id(0);
          const auto pointIndex = item.get_local_id(0);

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
                                                                             item);

          if (pointIndex == 0)
            devHasConverged[ltsFace] = hasConvergedLocal;

          devLocalSlipRate[ltsFace][pointIndex] =
              0.5 * (localSlipRateMagnitude + std::fabs(slipRateTest));
          devSlipRateMagnitude[ltsFace][pointIndex] = std::fabs(slipRateTest);

          devMu[ltsFace][pointIndex] = Derived::updateMu(
              localSlipRateMagnitude, localStateVariable, details, ltsFace, pointIndex);
        });
      });
    }
  }

  void calcSlipRateAndTraction(unsigned timeIndex) {
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

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

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
      });
    });
  }

  void saveDynamicStressOutput() {
    auto fullUpdateTime{this->mFullUpdateTime};
    auto muW{this->drParameters->muW};
    auto rsF0{this->drParameters->rsF0};

    auto* devDynStressTime{this->dynStressTime};
    auto* devDynStressTimePending{this->dynStressTimePending};
    auto* devRuptureTime{this->ruptureTime};
    auto* devMu{this->mu};

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        const auto localRuptureTime = devRuptureTime[ltsFace][pointIndex];
        if (localRuptureTime > 0.0 && localRuptureTime <= fullUpdateTime &&
            devDynStressTimePending[ltsFace][pointIndex] &&
            devMu[ltsFace][pointIndex] <= (muW + 0.05 * (rsF0 - muW))) {
          devDynStressTime[ltsFace][pointIndex] = fullUpdateTime;
          devDynStressTimePending[ltsFace][pointIndex] = false;
        }
      });
    });
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
                                      sycl::nd_item<1> item) {

    const auto ltsFace = item.get_group().get_group_id(0);
    const auto pointIndex = item.get_local_id(0);

    // Note that we need double precision here, since single precision led to NaNs.
    double muF{0.0}, dMuF{0.0}, g{0.0}, dG{0.0};
    slipRateTest = slipRateMagnitude;

    for (unsigned i = 0; i < solverSettings.maxNumberSlipRateUpdates; i++) {
      muF = Derived::updateMu(slipRateTest, localStateVariable, details, ltsFace, pointIndex);
      dMuF = Derived::updateMuDerivative(
          slipRateTest, localStateVariable, details, ltsFace, pointIndex);

      g = -invEtaS * (sycl::fabs(normalStress) * muF - absoluteShearStress) - slipRateTest;

      auto group = item.get_group();
      const bool converged =
          sycl::all_of_group(group, std::fabs(g) < solverSettings.newtonTolerance);

      if (converged)
        return true;

      dG = -invEtaS * (std::fabs(normalStress) * dMuF) - 1.0;
      slipRateTest =
          sycl::max(friction_law::rs::almostZero(), static_cast<real>(slipRateTest - (g / dG)));
    }
    return false;
  }

  void updateNormalStress(size_t timeIndex) {
    auto* devFaultStresses{this->faultStresses};
    auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
    auto* devNormalStress{this->initialVariables.normalStress};

    auto tpCurrentLayerDetails = tpMethod.getCurrentLayerDetails();

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);
        auto& faultStresses = devFaultStresses[ltsFace];

        devNormalStress[ltsFace][pointIndex] =
            std::min(static_cast<real>(0.0),
                     faultStresses.normalStress[timeIndex][pointIndex] +
                         devInitialStressInFaultCS[ltsFace][pointIndex][0] -
                         TPMethod::getFluidPressure(tpCurrentLayerDetails, ltsFace, pointIndex));
      });
    });
  }

  protected:
  real (*a)[misc::NumPaddedPoints];
  real (*sl0)[misc::NumPaddedPoints];
  real (*stateVariable)[misc::NumPaddedPoints];
  bool* hasConverged{};

  TPMethod tpMethod;
  rs::Settings settings{};
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_RATEANDSTATE_H
