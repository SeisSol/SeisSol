// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

namespace seissol::dr::friction_law::gpu {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived, class TPMethod>
class RateAndStateBase : public BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>> {
  public:
  explicit RateAndStateBase(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>>::BaseFrictionSolver(drParameters) {}

  ~RateAndStateBase() = default;

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {
    const auto* concreteLts = dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup);
    data->a = layerData.var(concreteLts->rsA, seissol::initializer::AllocationPlace::Device);
    data->sl0 = layerData.var(concreteLts->rsSl0, seissol::initializer::AllocationPlace::Device);
    data->stateVariable =
        layerData.var(concreteLts->stateVariable, seissol::initializer::AllocationPlace::Device);
    Derived::copySpecificLtsDataTreeToLocal(data, layerData, dynRup, fullUpdateTime);
    TPMethod::copyLtsTreeToLocal(data, layerData, dynRup, fullUpdateTime);
  }

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& ctx, unsigned timeIndex) {
    // compute initial slip rate and reference values
    Derived::calcInitialVariables(ctx, timeIndex);

    updateStateVariableIterative(ctx, timeIndex);
    Derived::executeIfNotConverged();

    TPMethod::calcFluidPressure(ctx, timeIndex, true);
    updateNormalStress(ctx, timeIndex);
    calcSlipRateAndTraction(ctx, timeIndex);
  }

  SEISSOL_DEVICE static void preHook(FrictionLawContext& ctx) {
    // copy state variable from last time step
    ctx.stateVariableBuffer = ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex];
  }

  SEISSOL_DEVICE static void postHook(FrictionLawContext& ctx) { Derived::resampleStateVar(ctx); }

  /*
   * Compute shear stress magnitude, localSlipRate, effective normal stress, reference state
   * variable. Also sets slipRateMagnitude member to reference value.
   */
  SEISSOL_DEVICE static void calcInitialVariables(FrictionLawContext& ctx, unsigned int timeIndex) {
    auto& devStateVariableBuffer{ctx.stateVariableBuffer};
    auto& devFaultStresses{ctx.faultStresses};
    auto* devSlipRateMagnitude{ctx.data->slipRateMagnitude};
    auto* devSlipRate1{ctx.data->slipRate1};
    auto* devSlipRate2{ctx.data->slipRate2};
    auto* devInitialStressInFaultCS{ctx.data->initialStressInFaultCS};

    auto& devAbsoluteShearTraction{ctx.initialVariables.absoluteShearTraction};
    auto& devLocalSlipRate{ctx.initialVariables.localSlipRate};
    auto& devStateVarReference{ctx.initialVariables.stateVarReference};

    updateNormalStress(ctx, timeIndex);

    auto& faultStresses = devFaultStresses;

    devStateVarReference = devStateVariableBuffer;

    const real totalTraction1 = devInitialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][3] +
                                faultStresses.traction1[timeIndex];

    const real totalTraction2 = devInitialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][5] +
                                faultStresses.traction2[timeIndex];

    devAbsoluteShearTraction = misc::magnitude(totalTraction1, totalTraction2);
    auto localSlipRateMagnitude = misc::magnitude(devSlipRate1[ctx.ltsFace][ctx.pointIndex],
                                                  devSlipRate2[ctx.ltsFace][ctx.pointIndex]);

    localSlipRateMagnitude = std::max(rs::almostZero(), localSlipRateMagnitude);
    devSlipRateMagnitude[ctx.ltsFace][ctx.pointIndex] = localSlipRateMagnitude;
    devLocalSlipRate = localSlipRateMagnitude;
  }

  SEISSOL_DEVICE static void updateStateVariableIterative(FrictionLawContext& ctx,
                                                          unsigned timeIndex) {
    auto& devLocalSlipRate{ctx.initialVariables.localSlipRate};
    auto& devStateVariableBuffer{ctx.stateVariableBuffer};
    auto& devNormalStress{ctx.initialVariables.normalStress};
    auto& devAbsoluteShearStress{ctx.initialVariables.absoluteShearTraction};
    auto* devMu{ctx.data->mu};

    rs::Settings settings{};

    for (unsigned j = 0; j < settings.numberStateVariableUpdates; j++) {

      const auto dt{ctx.data->deltaT[timeIndex]};
      Derived::updateStateVariable(ctx, dt);
      TPMethod::calcFluidPressure(ctx, timeIndex, false);
      updateNormalStress(ctx, timeIndex);

      const auto localStateVariable = devStateVariableBuffer;
      const auto normalStress = devNormalStress;
      const auto absoluteShearStress = devAbsoluteShearStress;
      const auto localSlipRateMagnitude = ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];
      const auto& localImpAndEta = ctx.data->impAndEta[ctx.ltsFace];

      auto& exportMu = devMu[ctx.ltsFace][ctx.pointIndex];

      real slipRateTest{};
      bool hasConvergedLocal = RateAndStateBase::invertSlipRateIterative(ctx,
                                                                         slipRateTest,
                                                                         localStateVariable,
                                                                         normalStress,
                                                                         absoluteShearStress,
                                                                         localSlipRateMagnitude,
                                                                         localImpAndEta.invEtaS,
                                                                         exportMu,
                                                                         settings);
      deviceBarrier(ctx);

      devLocalSlipRate = 0.5 * (localSlipRateMagnitude + std::fabs(slipRateTest));
      ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] = std::fabs(slipRateTest);
    }
  }

  SEISSOL_DEVICE static void calcSlipRateAndTraction(FrictionLawContext& ctx, unsigned timeIndex) {
    auto& devStateVarReference{ctx.initialVariables.stateVarReference};
    auto& devLocalSlipRate{ctx.initialVariables.localSlipRate};
    auto& devStateVariableBuffer{ctx.stateVariableBuffer}; // localStateVariable
    auto& devNormalStress{ctx.initialVariables.normalStress};
    auto& devAbsoluteTraction{ctx.initialVariables.absoluteShearTraction};
    auto& devFaultStresses{ctx.faultStresses};
    auto& devTractionResults{ctx.tractionResults};

    const auto deltaTime{ctx.data->deltaT[timeIndex]};

    Derived::updateStateVariable(ctx, ctx.data->deltaT[timeIndex]);

    const auto localStateVariable = devStateVariableBuffer;
    const auto slipRateMagnitude = ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];

    // the only mu calculation left, outside of the fixed-point loop
    auto details = Derived::getMuDetails(ctx, localStateVariable);
    ctx.data->mu[ctx.ltsFace][ctx.pointIndex] = Derived::updateMu(ctx, slipRateMagnitude, details);

    const real strength = -ctx.data->mu[ctx.ltsFace][ctx.pointIndex] * devNormalStress;

    const auto* initialStressInFaultCS =
        ctx.data->initialStressInFaultCS[ctx.ltsFace][ctx.pointIndex];
    const auto savedTraction1 = devFaultStresses.traction1[timeIndex];
    const auto savedTraction2 = devFaultStresses.traction2[timeIndex];

    // calculate absolute value of stress in Y and Z direction
    const real totalTraction1 = initialStressInFaultCS[3] + savedTraction1;
    const real totalTraction2 = initialStressInFaultCS[5] + savedTraction2;

    // Compute slip
    ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex] +=
        slipRateMagnitude * deltaTime;

    // Update slip rate
    const auto etaS = ctx.data->impAndEta[ctx.ltsFace].etaS;
    const auto divisor = strength + etaS * slipRateMagnitude;
    auto slipRate1 = slipRateMagnitude * totalTraction1 / divisor;
    auto slipRate2 = slipRateMagnitude * totalTraction2 / divisor;

    // calculate traction
    const auto traction1 = savedTraction1 - etaS * slipRate1;
    const auto traction2 = savedTraction2 - etaS * slipRate2;

    // Save traction for flux computation
    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] = traction1;
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] = traction2;

    // update directional slip
    ctx.data->slip1[ctx.ltsFace][ctx.pointIndex] += slipRate1 * deltaTime;
    ctx.data->slip2[ctx.ltsFace][ctx.pointIndex] += slipRate2 * deltaTime;

    // update traction
    devTractionResults.traction1[timeIndex] = traction1;
    devTractionResults.traction2[timeIndex] = traction2;

    // update slip rate
    ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] = slipRate1;
    ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] = slipRate2;
  }

  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& ctx) {
    auto fullUpdateTime{ctx.data->mFullUpdateTime};
    auto muW{ctx.data->drParameters.muW};
    auto rsF0{ctx.data->drParameters.rsF0};

    const auto localRuptureTime = ctx.data->ruptureTime[ctx.ltsFace][ctx.pointIndex];
    if (localRuptureTime > 0.0 && localRuptureTime <= fullUpdateTime &&
        ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] &&
        ctx.data->mu[ctx.ltsFace][ctx.pointIndex] <= (muW + 0.05 * (rsF0 - muW))) {
      ctx.data->dynStressTime[ctx.ltsFace][ctx.pointIndex] = fullUpdateTime;
      ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] = false;
    }
  }

  SEISSOL_DEVICE static bool invertSlipRateIterative(FrictionLawContext& ctx,
                                                     real& slipRateTest,
                                                     real localStateVariable,
                                                     real normalStress,
                                                     real absoluteShearStress,
                                                     real slipRateMagnitude,
                                                     real invEtaS,
                                                     real& exportMu,
                                                     rs::Settings solverSettings) {

    // Note that we need double precision here, since single precision led to NaNs.
    double muF{0.0};
    double dMuF{0.0};
    double g{0.0};
    double dG{0.0};
    slipRateTest = slipRateMagnitude;

    auto details = Derived::getMuDetails(ctx, localStateVariable);

    for (unsigned i = 0; i < solverSettings.maxNumberSlipRateUpdates; i++) {
      muF = Derived::updateMu(ctx, slipRateTest, details);

      g = -invEtaS * (std::fabs(normalStress) * muF - absoluteShearStress) - slipRateTest;

      const bool converged = std::fabs(g) < solverSettings.newtonTolerance;

      if (converged) {
        // we've reached the fixed point
        // NOTE: in doubt, a fixed-point mu can be recovered from slipRateTest at this point.
        // just invert -invEtaS * (std::fabs(normalStress) * muF - absoluteShearStress) ==
        // slipRateTest for muF in that case.
        exportMu = muF;
        return true;
      }

      dMuF = Derived::updateMuDerivative(ctx, slipRateTest, details);
      dG = -invEtaS * (std::fabs(normalStress) * dMuF) - 1.0;
      slipRateTest =
          std::max(friction_law::rs::almostZero(), static_cast<real>(slipRateTest - (g / dG)));
    }
    return false;
  }

  SEISSOL_DEVICE static void updateNormalStress(FrictionLawContext& ctx, size_t timeIndex) {
    ctx.initialVariables.normalStress =
        std::min(static_cast<real>(0.0),
                 ctx.faultStresses.normalStress[timeIndex] +
                     ctx.data->initialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][0] -
                     TPMethod::getFluidPressure(ctx));
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_
