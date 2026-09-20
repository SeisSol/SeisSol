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

#include <cmath>
#include <limits>

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

    constexpr auto Place = seissol::initializer::AllocationPlace::Device;

    data->a = layerData.var<LTSRateAndState::RsA>(Place);
    data->sl0 = layerData.var<LTSRateAndState::RsSl0>(Place);
    data->stateVariable = layerData.var<LTSRateAndState::StateVariable>(Place);
    data->f0 = layerData.var<LTSRateAndState::RsF0>(Place);
    data->muW = layerData.var<LTSRateAndState::RsMuW>(Place);
    data->b = layerData.var<LTSRateAndState::RsB>(Place);
    data->convergenceInner = layerData.var<LTSRateAndState::ConvergenceInner>(Place);
    data->convergenceOuter = layerData.var<LTSRateAndState::ConvergenceOuter>(Place);

    Derived::copySpecificStorageDataToLocal(data, layerData);
    TPMethod::copyStorageToLocal(data, layerData);
  }

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& __restrict ctx,
                                                   uint32_t timeIndex) {
    // compute initial slip rate and reference values
    Derived::calcInitialVariables(ctx);

    updateStateVariableIterative(ctx, timeIndex);

    TPMethod::calcFluidPressure(ctx, timeIndex, true);
    updateDirectionAndProjections(ctx);
    updateNormalStress(ctx);
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
  SEISSOL_DEVICE static void calcInitialVariables(FrictionLawContext& __restrict ctx) {
    ctx.initialVariables.stateVarReference = ctx.stateVariableBuffer;

    const real totalTraction1 = ctx.data->initialStressInFaultCS[ctx.ltsFace][3][ctx.pointIndex] +
                                ctx.faultStresses.traction1;

    const real totalTraction2 = ctx.data->initialStressInFaultCS[ctx.ltsFace][5][ctx.pointIndex] +
                                ctx.faultStresses.traction2;

    ctx.initialVariables.absoluteShearTraction = misc::magnitude(totalTraction1, totalTraction2);

    ctx.initialVariables.etaNormal =
        common::projectEtaNormal(ctx.data->impAndEta[ctx.ltsFace],
                                 ctx.data->impedanceMatrices[ctx.ltsFace],
                                 totalTraction1,
                                 totalTraction2,
                                 ctx.initialVariables.absoluteShearTraction);

    // initial slip direction: the trial traction. For isotropy this stays exact.
    const real invAbsolute =
        (ctx.initialVariables.absoluteShearTraction > 0)
            ? static_cast<real>(1.0) / ctx.initialVariables.absoluteShearTraction
            : static_cast<real>(0.0);
    ctx.initialVariables.slipDirection1 = totalTraction1 * invAbsolute;
    ctx.initialVariables.slipDirection2 = totalTraction2 * invAbsolute;

    auto localSlipRateMagnitude = misc::magnitude(ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex],
                                                  ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex]);

    localSlipRateMagnitude = std::max(rs::almostZero(), localSlipRateMagnitude);
    ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] = localSlipRateMagnitude;
    ctx.initialVariables.localSlipRate = localSlipRateMagnitude;

    // after the slip rate is set: updateNormalStress reads it (and etaNormal)
    updateNormalStress(ctx);
  }

  /**
   * Anisotropy: sweeps the slip direction, and with it every direction dependent projection of the
   * impedance. A no-op for isotropy, where the slip is exactly parallel to the trial traction.
   *
   * The exact condition is tau0 = (S I + V eta_ss) n with |n| = 1, see common::updateSlipDirection.
   * The strength belonging to the current slip rate is recovered as S = n^T tau0 - V * eta_proj,
   * so no additional friction law evaluation is needed. Riding along with the existing outer
   * fixed-point loop makes the sweep essentially free.
   *
   * @returns 1 / eta_proj for the (possibly updated) direction
   */
  SEISSOL_DEVICE static real updateDirectionAndProjections(FrictionLawContext& __restrict ctx) {
    const real totalTraction1 = ctx.data->initialStressInFaultCS[ctx.ltsFace][3][ctx.pointIndex] +
                                ctx.faultStresses.traction1;

    const real totalTraction2 = ctx.data->initialStressInFaultCS[ctx.ltsFace][5][ctx.pointIndex] +
                                ctx.faultStresses.traction2;

    if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
      const auto [etaProj, unusedInv] = common::projectEta(ctx.data->impAndEta[ctx.ltsFace],
                                                           ctx.data->impedanceMatrices[ctx.ltsFace],
                                                           ctx.initialVariables.slipDirection1,
                                                           ctx.initialVariables.slipDirection2,
                                                           static_cast<real>(1.0));

      const real slipRate = ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];
      const real strength = ctx.initialVariables.absoluteShearTraction - slipRate * etaProj;

      const auto [n1, n2] =
          common::updateSlipDirection(ctx.data->impAndEta[ctx.ltsFace],
                                      ctx.data->impedanceMatrices[ctx.ltsFace],
                                      strength,
                                      slipRate,
                                      totalTraction1,
                                      totalTraction2,
                                      misc::magnitude(totalTraction1, totalTraction2));

      ctx.initialVariables.slipDirection1 = n1;
      ctx.initialVariables.slipDirection2 = n2;
      ctx.initialVariables.absoluteShearTraction = n1 * totalTraction1 + n2 * totalTraction2;
      ctx.initialVariables.etaNormal =
          common::projectEtaNormal(ctx.data->impAndEta[ctx.ltsFace],
                                   ctx.data->impedanceMatrices[ctx.ltsFace],
                                   n1,
                                   n2,
                                   static_cast<real>(1.0));

      const auto [unusedNewEta, newInvEta] =
          common::projectEta(ctx.data->impAndEta[ctx.ltsFace],
                             ctx.data->impedanceMatrices[ctx.ltsFace],
                             n1,
                             n2,
                             static_cast<real>(1.0));
      return newInvEta;
    } else {
      const auto [unusedEta, invEta] =
          common::projectEta(ctx.data->impAndEta[ctx.ltsFace],
                             ctx.data->impedanceMatrices[ctx.ltsFace],
                             totalTraction1,
                             totalTraction2,
                             ctx.initialVariables.absoluteShearTraction);
      return invEta;
    }
  }

  SEISSOL_DEVICE static void updateStateVariableIterative(FrictionLawContext& __restrict ctx,
                                                          uint32_t timeIndex) {
    bool hasConvergedOuter = false;
    bool hasConvergedInner = true;

    for (uint32_t j = 0; j < ctx.data->drParameters.rsNumberStateVariableUpdates; j++) {

      const auto dt{ctx.args->deltaT[timeIndex]};
      Derived::updateStateVariable(ctx, dt);
      TPMethod::calcFluidPressure(ctx, timeIndex, false);
      const real invEta = updateDirectionAndProjections(ctx);
      updateNormalStress(ctx);

      const auto localStateVariable = ctx.stateVariableBuffer;
      const auto normalStress = ctx.initialVariables.normalStress;
      const auto normalStressStick = ctx.initialVariables.normalStressStick;
      const auto absoluteShearStress = ctx.initialVariables.absoluteShearTraction;
      const auto localSlipRateMagnitude = ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];

      real slipRateTest{0};
      real exportMu{0};

      const bool hasConvergedLocal =
          RateAndStateBase::invertSlipRateIterative(ctx,
                                                    slipRateTest,
                                                    localStateVariable,
                                                    normalStress,
                                                    normalStressStick,
                                                    ctx.initialVariables.etaNormal,
                                                    absoluteShearStress,
                                                    localSlipRateMagnitude,
                                                    invEta,
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

    const auto savedTraction1 = ctx.faultStresses.traction1;
    const auto savedTraction2 = ctx.faultStresses.traction2;

    // Compute slip
    ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex] +=
        slipRateMagnitude * deltaTime;

    // the direction along which the slip rate is decomposed; scaled such that dividing by
    // `divisor` yields the unit slip direction. For isotropy slipDirection is the normalized trial
    // traction and absoluteShearTraction its magnitude, so this is the previous expression.
    const real dirTraction1 =
        ctx.initialVariables.slipDirection1 * ctx.initialVariables.absoluteShearTraction;
    const real dirTraction2 =
        ctx.initialVariables.slipDirection2 * ctx.initialVariables.absoluteShearTraction;

    const auto [etaS, _] = common::projectEta(ctx.data->impAndEta[ctx.ltsFace],
                                              ctx.data->impedanceMatrices[ctx.ltsFace],
                                              ctx.initialVariables.slipDirection1,
                                              ctx.initialVariables.slipDirection2,
                                              static_cast<real>(1.0));

    // Update slip rate
    const auto divisor = strength + etaS * slipRateMagnitude;
    const auto slipRate1 = slipRateMagnitude * dirTraction1 / divisor;
    const auto slipRate2 = slipRateMagnitude * dirTraction2 / divisor;

    const auto [tU1, tU2] = common::matmulEta(ctx.data->impAndEta[ctx.ltsFace],
                                              ctx.data->impedanceMatrices[ctx.ltsFace],
                                              slipRate1,
                                              slipRate2);

    // calculate traction
    const auto traction1 = savedTraction1 - tU1;
    const auto traction2 = savedTraction2 - tU2;

    // Save traction for flux computation
    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] = traction1;
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] = traction2;

    // update directional slip
    ctx.data->slip1[ctx.ltsFace][ctx.pointIndex] += slipRate1 * deltaTime;
    ctx.data->slip2[ctx.ltsFace][ctx.pointIndex] += slipRate2 * deltaTime;

    // update traction
    // note that the normal stress written here is the *dynamic* normal traction, i.e. in the same
    // space as faultStresses/qInterpolated -- not the effective normal stress used for the friction
    // strength above (which additionally carries the initial stress, the initial pressure, thermal
    // pressurization and the min(., 0) clamp).
    ctx.tractionResults.normalStress =
        ctx.faultStresses.normalStress - slipRateMagnitude * ctx.initialVariables.etaNormal;
    ctx.tractionResults.traction1 = traction1;
    ctx.tractionResults.traction2 = traction2;

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

  /**
   * The effective normal stress a trial slip rate belongs to.
   *
   * Without the anisotropic normal coupling this is the value updateNormalStress left behind and
   * the solve is the one every other material runs. With it, sigma follows the slip rate, and
   * evaluating it inside the Newton moves the coupling out of the outer fixed point, which
   * resolves it at a linear rate, into the quadratic one.
   */
  SEISSOL_DEVICE static real effectiveNormalStress(real normalStress,
                                                   real normalStressStick,
                                                   real etaNormal,
                                                   real slipRate) {
    if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
      return std::min(static_cast<real>(0.0), normalStressStick - slipRate * etaNormal);
    } else {
      return normalStress;
    }
  }

  SEISSOL_DEVICE static bool invertSlipRateIterative(FrictionLawContext& __restrict ctx,
                                                     real& slipRateTest,
                                                     real localStateVariable,
                                                     real normalStress,
                                                     real normalStressStick,
                                                     real etaNormal,
                                                     real absoluteShearStress,
                                                     real slipRateMagnitude,
                                                     real invEtaS,
                                                     real& exportMu) {
    // Solve  g(V) = -invEtaS * (|sigma(V)| * mu(V) - tau) - V = 0   for V = slipRateTest,
    // with sigma following the slip rate through the anisotropic normal coupling. The bracket is
    // closed-form and survives that coupling: it only needs mu(0) = 0, mu >= 0, |sigma| >= 0
    // (no search, endpoints not evaluated):
    //   g(0+)          = invEtaS * tau             > 0
    //   g(tau*invEtaS) = -invEtaS*|sigma(.)|*mu(.) <= 0
    // Without the coupling dG < -1 everywhere and the root is unique; the coupling can weaken dG
    // (cf. below), and bisection converges to a root in the bracket either way.
    // rtsafe: Newton while it stays in the bracket and outruns bisection, else bisect.
    // The bracket is non-increasing and loses half of its decades on every fallback => the
    // iterate settles and termination is relative in V-space (|dV| < xacc * V), with two floors
    // that keep the test reachable in finite precision: xacc clamped to a few ulp, and a residual
    // that has sunk into the rounding noise of its own evaluation.

    const auto details = Derived::getMuDetails(ctx, localStateVariable);
    const real tau = absoluteShearStress;

    real xLow = friction_law::rs::almostZero();
    real xHigh = std::max(xLow, tau * invEtaS); // tau~0 => collapses to ~0, root ~0
    real x = std::min(std::max(slipRateMagnitude, xLow), xHigh); // warm start, clamped
    real dx = xHigh - xLow;                                      // becomes dxOld on first iteration

    // Number of roundings that enter one residual evaluation; used to size both floors below.
    constexpr real NoiseFactor = 4;
    constexpr real Eps = std::numeric_limits<real>::epsilon();

    // rsSlipRateTolerance is a RELATIVE step tolerance. A nonzero step is at least one ulp of the
    // iterate, so a tolerance below Eps cannot be met at all and would leave the exact-fixed-point
    // guard as the only way out. Clamp it to a few ulp.
    const real xacc =
        std::max(static_cast<real>(ctx.data->drParameters.rsSlipRateTolerance), NoiseFactor * Eps);

    real muF{0};
    bool converged = false;

    // A point that carries no normal stress at the free-slip limit has its root exactly there:
    // |sigma| vanishes, g is the line tau * invEtaS - V, and g(xHigh) = 0. That is worth taking
    // directly, because it is the one root rtsafe cannot approach: a root on the bracket boundary
    // leaves the Newton step the same size as the previous one, so the guard falls back to
    // bisection on every iteration and the solve spends its whole budget halving.
    if (effectiveNormalStress(normalStress, normalStressStick, etaNormal, xHigh) ==
        static_cast<real>(0.0)) {
      x = xHigh;
      converged = true;
    }

    for (uint32_t i = 0; i < ctx.data->drParameters.rsMaxNumberSlipRateUpdates; i++) {
      const bool active = !converged;

      // >>> precision knob: evaluate muF/g/dG in double (promote sigma, tau, x) to drop the
      //     noise floor AND make the sign below exact. Needs a double mu() evaluation.
      muF = Derived::updateMu(ctx, x, details);
      const real dMuF = Derived::updateMuDerivative(ctx, x, details);
      // sigma follows the trial slip rate, so it is evaluated at x rather than taken frozen: that
      // moves the normal coupling out of the outer fixed point and into this Newton.
      const real sigma = effectiveNormalStress(normalStress, normalStressStick, etaNormal, x);
      const real absSigma = std::abs(sigma);
      const real g = -invEtaS * (absSigma * muF - tau) - x;

      // |sigma| = -sigma while the fault is closed, and sigma follows the slip rate through the
      // anisotropic normal coupling, so d|sigma|/dV = etaNormal there.
      real dAbsSigma{};
      if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
        dAbsSigma = (sigma < static_cast<real>(0)) ? etaNormal : static_cast<real>(0);
      } else {
        dAbsSigma = static_cast<real>(0);
      }
      const real dGFrozen = -invEtaS * (absSigma * dMuF) - static_cast<real>(1);
      const real dGCoupled = -invEtaS * (absSigma * dMuF + dAbsSigma * muF) - static_cast<real>(1);
      // A fault that loses normal stress as it slips (etaNormal < 0) is the only case in which the
      // coupling can weaken g. It stays strictly decreasing as long as
      // |etaNormal| * mu < eta_proj + |sigma| * mu', which is NOT a comfortable margin: a locked
      // point at |etaNormal| = eta_proj already peaks at dGCoupled = -0.06 over the bracket and
      // turns positive slightly above that ratio, a slipping one holds out to roughly three times
      // eta_proj. Positive definiteness bounds the ratio by sqrt(eta_nn / eta_proj), so a strongly
      // anisotropic material lands close to the edge. Past it g has several roots in the bracket
      // and the solve returns one of them, which is why the bracket rather than the derivative
      // carries the robustness. dGFrozen is negative by construction, which keeps the Newton step
      // pointing into the bracket.
      const real dG = (dGCoupled < static_cast<real>(0)) ? dGCoupled : dGFrozen;

      // |sigma| * mu and tau cancel at the root, so the rounding error of g does not shrink with
      // the iterate: it stays at Eps times the magnitude of the two cancelling terms. Below that
      // level the sign of g -- and with it the bracket update -- carries no information.
      const real gNoise = NoiseFactor * Eps * invEtaS * (absSigma * muF + tau);

      // maintain the straddling bracket from sign(g) (g decreasing):
      //   g > 0 => root at larger  V => raise lower bound
      //   g < 0 => root at smaller V => lower upper bound
      const bool gPos = g > static_cast<real>(0);
      xLow = (active && gPos) ? x : xLow;
      xHigh = (active && !gPos) ? x : xHigh;

      const real dxOld = dx;
      const real xNewton = x - g / dG;
      // Bisect geometrically. The bracket spans the whole admissible range of slip rates, from
      // almostZero() up to the free-slip limit tau/eta_s, so its arithmetic midpoint sits many
      // orders of magnitude above the root of a locked or creeping point, and a fallback would
      // then need one halving per factor of two to walk back down. The geometric midpoint halves
      // the number of decades instead, which is the scale the root lives on. The two square roots
      // keep the product from underflowing for the smallest brackets.
      const real xBisect = std::sqrt(xLow) * std::sqrt(xHigh);

      // bisect if Newton leaves the bracket or does not outrun bisection
      const bool useBisect = (xNewton <= xLow) || (xNewton >= xHigh) ||
                             (std::abs(static_cast<real>(2) * g) > std::abs(dxOld * dG));
      const real xUpdated = useBisect ? xBisect : xNewton;
      dx = xUpdated - x;

      // V-space convergence, residual-noise floor, and the no-representable-change guard
      // => cannot livelock
      const bool nowConverged =
          (std::abs(dx) < xacc * std::abs(x)) || (std::abs(g) <= gNoise) || (xUpdated == x);

      // advance active, not-yet-converged lanes; freeze at the converged point so that
      // slipRateTest and exportMu are the mu-consistent pair at that point
      converged |= nowConverged;
      x = converged ? x : xUpdated;

      deviceWarpBarrier(ctx);
      if (deviceWarpAll(ctx, converged)) {
        break;
      }
    }

    slipRateTest = x;
    exportMu = Derived::updateMu(ctx, x, details);
    return converged;
  }

  /**
   * Effective normal stress, including the anisotropic normal/shear coupling.
   *
   * With the 3x3 impedance eta, shear slip changes the fault-normal traction:
   *   sigma(V) = sigma_stick - V * etaNormal,   etaNormal = (eta * n)_n
   * For every isotropic material etaNormal is zero and this reduces to the previous formula.
   *
   * The slip rate is taken from the previous outer fixed-point iteration (or, on entry, from the
   * previous time step). normalStressStick keeps the part that does not depend on it, so that the
   * Newton solve can follow sigma(V) itself.
   */
  SEISSOL_DEVICE static void updateNormalStress(FrictionLawContext& __restrict ctx) {
    ctx.initialVariables.normalStressStick =
        ctx.faultStresses.normalStress +
        ctx.data->initialStressInFaultCS[ctx.ltsFace][0][ctx.pointIndex] +
        ctx.faultStresses.fluidPressure + ctx.data->initialPressure[ctx.ltsFace][ctx.pointIndex] -
        TPMethod::getFluidPressure(ctx);
    ctx.initialVariables.normalStress =
        std::min(static_cast<real>(0.0),
                 ctx.initialVariables.normalStressStick -
                     ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] *
                         ctx.initialVariables.etaNormal);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_
