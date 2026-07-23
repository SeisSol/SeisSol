// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::friction_law::gpu {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived, class TPMethod>
class RateAndStateBase : public BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>> {
  public:
  explicit RateAndStateBase(const FrictionLawParameters& drParameters)
      : BaseFrictionSolver<RateAndStateBase<Derived, TPMethod>>::BaseFrictionSolver(drParameters) {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  ~RateAndStateBase() override = default;

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {
    data->a = layerData.var<LTSRateAndState::RsA>(seissol::initializer::AllocationPlace::Device);
    data->sl0 =
        layerData.var<LTSRateAndState::RsSl0>(seissol::initializer::AllocationPlace::Device);
    data->stateVariable = layerData.var<LTSRateAndState::StateVariable>(
        seissol::initializer::AllocationPlace::Device);
    data->f0 = layerData.var<LTSRateAndState::RsF0>(seissol::initializer::AllocationPlace::Device);
    data->muW =
        layerData.var<LTSRateAndState::RsMuW>(seissol::initializer::AllocationPlace::Device);
    data->b = layerData.var<LTSRateAndState::RsB>(seissol::initializer::AllocationPlace::Device);
    data->convergenceInner = layerData.var<LTSRateAndState::ConvergenceInner>();
    data->convergenceOuter = layerData.var<LTSRateAndState::ConvergenceOuter>();
    Derived::copySpecificStorageDataToLocal(data, layerData);
    TPMethod::copyStorageToLocal(data, layerData);
  }

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& __restrict ctx,
                                                   uint32_t timeIndex) {
    // compute initial slip rate and reference values
    Derived::calcInitialVariables(ctx, timeIndex);

    updateStateVariableIterative(ctx, timeIndex);

    TPMethod::calcFluidPressure(ctx, timeIndex, true);
    updateNormalStress(ctx, timeIndex);
    calcSlipRateAndTraction(ctx, timeIndex);
  }

  SEISSOL_DEVICE static void preHook(FrictionLawContext& __restrict ctx) {
    // copy state variable from last time step
    ctx.stateVariableBuffer = ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex];
  }

  SEISSOL_DEVICE static void postHook(FrictionLawContext& __restrict ctx) {
    Derived::resampleStateVar(ctx);
  }

  /*
   * Compute shear stress magnitude, localSlipRate, effective normal stress, reference state
   * variable. Also sets slipRateMagnitude member to reference value.
   */
  SEISSOL_DEVICE static void calcInitialVariables(FrictionLawContext& __restrict ctx,
                                                  uint32_t timeIndex) {
    updateNormalStress(ctx, timeIndex);

    ctx.initialVariables.stateVarReference = ctx.stateVariableBuffer;

    const real totalTraction1 = ctx.data->initialStressInFaultCS[ctx.ltsFace][3][ctx.pointIndex] +
                                ctx.faultStresses.traction1[timeIndex];

    const real totalTraction2 = ctx.data->initialStressInFaultCS[ctx.ltsFace][5][ctx.pointIndex] +
                                ctx.faultStresses.traction2[timeIndex];

    ctx.initialVariables.absoluteShearTraction = misc::magnitude(totalTraction1, totalTraction2);
    auto localSlipRateMagnitude = misc::magnitude(ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex],
                                                  ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex]);

    localSlipRateMagnitude = std::max(rs::almostZero(), localSlipRateMagnitude);
    ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] = localSlipRateMagnitude;
    ctx.initialVariables.localSlipRate = localSlipRateMagnitude;
  }

  SEISSOL_DEVICE static void updateStateVariableIterative(FrictionLawContext& __restrict ctx,
                                                          uint32_t timeIndex) {

    bool hasConvergedOuter = false;
    bool hasConvergedInner = true;

    for (uint32_t j = 0; j < ctx.data->drParameters.rsNumberStateVariableUpdates; j++) {

      const auto dt{ctx.args->deltaT[timeIndex]};
      Derived::updateStateVariable(ctx, dt);
      TPMethod::calcFluidPressure(ctx, timeIndex, false);
      updateNormalStress(ctx, timeIndex);

      const auto localStateVariable = ctx.stateVariableBuffer;
      const auto normalStress = ctx.initialVariables.normalStress;
      const auto absoluteShearStress = ctx.initialVariables.absoluteShearTraction;
      const auto localSlipRateMagnitude = ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];
      const auto& localImpAndEta = ctx.data->impAndEta[ctx.ltsFace];

      real slipRateTest{0};
      real exportMu{0};

      const bool hasConvergedLocal =
          RateAndStateBase::invertSlipRateIterative(ctx,
                                                    slipRateTest,
                                                    localStateVariable,
                                                    normalStress,
                                                    absoluteShearStress,
                                                    localSlipRateMagnitude,
                                                    localImpAndEta.invEtaS,
                                                    exportMu);

      hasConvergedInner &= hasConvergedLocal;

      ctx.initialVariables.localSlipRate = (localSlipRateMagnitude + slipRateTest) / 2;
      ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] = slipRateTest;
      ctx.data->mu[ctx.ltsFace][ctx.pointIndex] = exportMu;

      hasConvergedOuter =
          std::abs(localSlipRateMagnitude - slipRateTest) < ctx.data->drParameters.rsStateTolerance;

      // exit early and prevent thread/load data divergence
      deviceWarpBarrier(ctx);
      if (deviceWarpAll(ctx, hasConvergedOuter)) {
        break;
      }
    }
    deviceBarrier(ctx);
    ctx.data->convergenceOuter[ctx.ltsFace][ctx.pointIndex] &= hasConvergedOuter;
    ctx.data->convergenceInner[ctx.ltsFace][ctx.pointIndex] &= hasConvergedInner;
  }

  SEISSOL_DEVICE static void calcSlipRateAndTraction(FrictionLawContext& __restrict ctx,
                                                     uint32_t timeIndex) {
    const auto deltaTime{ctx.args->deltaT[timeIndex]};

    Derived::updateStateVariable(ctx, deltaTime);

    const auto localStateVariable = ctx.stateVariableBuffer;
    const auto slipRateMagnitude = ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];

    // the only mu calculation left, outside of the fixed-point loop
    const auto details = Derived::getMuDetails(ctx, localStateVariable);
    const auto mu = Derived::updateMu(ctx, slipRateMagnitude, details);

    ctx.data->mu[ctx.ltsFace][ctx.pointIndex] = mu;

    const real strength = -mu * ctx.initialVariables.normalStress;

    const auto* initialStressInFaultCS = ctx.data->initialStressInFaultCS[ctx.ltsFace];
    const auto savedTraction1 = ctx.faultStresses.traction1[timeIndex];
    const auto savedTraction2 = ctx.faultStresses.traction2[timeIndex];

    // calculate absolute value of stress in Y and Z direction
    const real totalTraction1 = initialStressInFaultCS[3][ctx.pointIndex] + savedTraction1;
    const real totalTraction2 = initialStressInFaultCS[5][ctx.pointIndex] + savedTraction2;

    // Compute slip
    ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex] +=
        slipRateMagnitude * deltaTime;

    // Update slip rate
    const auto etaS = ctx.data->impAndEta[ctx.ltsFace].etaS;
    const auto divisor = strength + etaS * slipRateMagnitude;
    const auto slipRate1 = slipRateMagnitude * totalTraction1 / divisor;
    const auto slipRate2 = slipRateMagnitude * totalTraction2 / divisor;

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
    ctx.tractionResults.traction1[timeIndex] = traction1;
    ctx.tractionResults.traction2[timeIndex] = traction2;

    // update slip rate
    ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] = slipRate1;
    ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] = slipRate2;
  }

  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& __restrict ctx,
                                                     real time) {
    auto muW{ctx.data->muW[ctx.ltsFace][ctx.pointIndex]};
    auto rsF0{ctx.data->f0[ctx.ltsFace][ctx.pointIndex]};

    const auto localRuptureTime = ctx.data->ruptureTime[ctx.ltsFace][ctx.pointIndex];
    if (localRuptureTime > static_cast<real>(0.0) && localRuptureTime <= time &&
        ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] &&
        ctx.data->mu[ctx.ltsFace][ctx.pointIndex] <=
            (muW + static_cast<real>(0.05) * (rsF0 - muW))) {
      ctx.data->dynStressTime[ctx.ltsFace][ctx.pointIndex] = time;
      ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] = false;
    }
  }

  SEISSOL_DEVICE static bool invertSlipRateIterative(FrictionLawContext& __restrict ctx,
                                                     real& slipRateTest,
                                                     real localStateVariable,
                                                     real normalStress,
                                                     real absoluteShearStress,
                                                     real slipRateMagnitude,
                                                     real invEtaS,
                                                     real& exportMu) {
    // Solve  g(V) = -invEtaS * (|sigma_n| * mu(V) - tau) - V = 0   for V = slipRateTest.
    // dG = -invEtaS*|sigma_n|*mu'(V) - 1 < -1 everywhere (mu' > 0), so g is strictly
    // decreasing: unique root, closed-form bracket (no search, endpoints not evaluated):
    //   g(0+)          = invEtaS * tau            > 0
    //   g(tau*invEtaS) = -invEtaS*|sigma_n|*mu(.) < 0
    // rtsafe: Newton while it stays in the bracket and outruns bisection, else bisect.
    // Bracket width is non-increasing and halves on every fallback => termination is in
    // V-space (|dV| < xacc), immune to the residual noise floor that makes a |g|-threshold
    // circulate at low precision.

    const auto details = Derived::getMuDetails(ctx, localStateVariable);
    const real absN = std::abs(normalStress);
    const real tau = absoluteShearStress;

    real xLow = friction_law::rs::almostZero();
    real xHigh = std::max(xLow, tau * invEtaS); // tau~0 => collapses to ~0, root ~0
    real x = std::min(std::max(slipRateMagnitude, xLow), xHigh); // warm start, clamped
    real dx = xHigh - xLow;                                      // becomes dxOld on first iteration

    // rsSlipRateTolerance is dimensionally a slip rate => reuse directly as the V-space
    // step tolerance. Reproduces the current convergence scale, but now terminating.
    const real xacc = ctx.data->drParameters.rsSlipRateTolerance;

    real muF{0};
    bool converged = false;

    for (uint32_t i = 0; i < ctx.data->drParameters.rsMaxNumberSlipRateUpdates; i++) {
      const bool active = !converged;

      // >>> precision knob: evaluate muF/g/dG in double (promote absN, tau, x) to drop the
      //     noise floor AND make the sign below exact. Needs a double mu() evaluation.
      muF = Derived::updateMu(ctx, x, details);
      const real dMuF = Derived::updateMuDerivative(ctx, x, details);
      const real g = -invEtaS * (absN * muF - tau) - x;
      const real dG = -invEtaS * (absN * dMuF) - static_cast<real>(1);

      // maintain the straddling bracket from sign(g) (g decreasing):
      //   g > 0 => root at larger  V => raise lower bound
      //   g < 0 => root at smaller V => lower upper bound
      const bool gPos = g > static_cast<real>(0);
      xLow = (active && gPos) ? x : xLow;
      xHigh = (active && !gPos) ? x : xHigh;

      const real dxOld = dx;
      const real xNewton = x - g / dG;
      const real xBisect = static_cast<real>(0.5) * (xLow + xHigh);

      // bisect if Newton leaves the bracket or does not outrun bisection
      const bool useBisect = (xNewton <= xLow) || (xNewton >= xHigh) ||
                             (std::abs(static_cast<real>(2) * g) > std::abs(dxOld * dG));
      const real xUpdated = useBisect ? xBisect : xNewton;
      dx = xUpdated - x;

      // V-space convergence + no-representable-change guard => cannot livelock
      const bool nowConverged = (std::abs(dx) < xacc * std::abs(x)) || (xUpdated == x);

      // advance active, not-yet-converged lanes; freeze at the converged point so that
      // slipRateTest and exportMu are the mu-consistent pair at that point
      x = (active && !nowConverged) ? xUpdated : x;
      converged |= nowConverged;

      deviceWarpBarrier(ctx);
      if (deviceWarpAll(ctx, converged)) {
        break;
      }
    }

    slipRateTest = x;
    exportMu = muF;
    return converged;
  }

  SEISSOL_DEVICE static void updateNormalStress(FrictionLawContext& __restrict ctx,
                                                uint32_t timeIndex) {
    ctx.initialVariables.normalStress =
        std::min(static_cast<real>(0.0),
                 ctx.faultStresses.normalStress[timeIndex] +
                     ctx.data->initialStressInFaultCS[ctx.ltsFace][0][ctx.pointIndex] +
                     ctx.faultStresses.fluidPressure[timeIndex] +
                     ctx.data->initialPressure[ctx.ltsFace][ctx.pointIndex] -
                     TPMethod::getFluidPressure(ctx));
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_
