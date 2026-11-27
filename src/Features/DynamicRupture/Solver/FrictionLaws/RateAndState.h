// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_RATEANDSTATE_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_RATEANDSTATE_H_

#include "Core/Common/Marker.h"
#include "Core/Common/MaterialType.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Physics/Dual.h"
#include "Features/DynamicRupture/Physics/RateAndStateCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#ifdef __INTEL_LLVM_COMPILER
#if __INTEL_LLVM_COMPILER >= 20250000
#define SEISSOL_INTEL_SIMD_EXCEPTION
// icpx 2026.0 still crashes in its vectorizer (vplan-vec) in single precision, on the strength
// and slip rate loop of calcSlipRateAndTraction and on the loops of invertSlipRateIterative; they
// stay scalar for it in either precision. Re-check with newer releases.
#if __INTEL_LLVM_COMPILER < 20260100
#define SEISSOL_INTEL_SIMD_EXCEPTION_STRICT
#endif
#endif
#endif // __INTEL_LLVM_COMPILER

namespace seissol::dr::friction_law {

namespace rs_detail {
/// whether the loops that some releases of icpx fail to vectorize (see above) are vectorized
#ifdef SEISSOL_INTEL_SIMD_EXCEPTION_STRICT
constexpr bool VectorizeStrict = false;
#else
constexpr bool VectorizeStrict = true;
#endif
#ifdef SEISSOL_INTEL_SIMD_EXCEPTION
constexpr bool VectorizeSelect = false;
#else
constexpr bool VectorizeSelect = true;
#endif
} // namespace rs_detail

/**
 * What a rate-and-state law keeps of the points of the range during a sub-step: the variables that
 * are computed initially in each sub-step and those that the fixed point of the slip rate sweeps.
 */
template <typename Cfg, common::RangeType Type>
struct RateAndStateStep {
  /// the state variable at the start of the step, which the state is integrated from. Careful, the
  /// state variable must always be corrected from it and not from the current state variable!
  PointValues<Cfg, Type> stateVarReference{};
  /// the shear traction along the slip direction; direction dependent, see
  /// updateDirectionAndProjections
  PointValues<Cfg, Type> absoluteShearTraction{};
  /// the slip rate the state variable is integrated with: the mean of the one at the start of the
  /// step and the latest iterate (Kaneko 2008, step 6)
  PointValues<Cfg, Type> localSlipRate{};
  PointValues<Cfg, Type> normalStress{};
  /// the same, before the slip rate dependent part and the clamp; the Newton solve needs it to
  /// evaluate sigma(V) itself
  PointValues<Cfg, Type> normalStressStick{};
  /// direction dependent, as the three below
  PointValues<Cfg, Type> etaInv{};
  PointValues<Cfg, Type> etaNormal{};
  /// unit slip direction; equals the normalized trial traction unless the impedance is anisotropic
  PointValues<Cfg, Type> slipDirection1{};
  PointValues<Cfg, Type> slipDirection2{};
};

/**
 * General implementation of a rate and state solver, once for the host and the device (see
 * FrictionLawDriver). The law `LawT` supplies, at a point of a face:
 *   StateScalar, FoldsStateVariable, MuDetails<Type>,
 *   copyLawStorageToLocal(data, layerData, place),
 *   stateVariableAt(ctx, pointIndex, stateVarReference, timeIncrement, slipRate),
 *   updateMuFolded(ctx, pointIndex, slipRate, stateVarReference, timeIncrement),
 *   getMuDetails(ctx, pointIndex, stateVariable, details), updateMu(ctx, pointIndex, slipRate,
 *   details), resampleStateVar(ctx).
 * The thermal pressurization `TPMethod` (NoTP, ThermalPressurization) supplies the fluid pressure.
 */
template <typename Cfg, typename LawT, typename TPMethod>
class RateAndState {
  private:
  // only the CRTP subclass `LawT` constructs this base
  RateAndState() = default;
  friend LawT;

  public:
  using real = Real<Cfg>;

  // the slip rate of a sub-step is the one the friction solve finds at its time point
  static constexpr bool SlipRateAtIntervalEnds = false;

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData,
                                             seissol::initializer::AllocationPlace place) {
    data->a = layerData.var<LTSRateAndState::RsA>(Cfg(), place);
    data->sl0 = layerData.var<LTSRateAndState::RsSl0>(Cfg(), place);
    data->stateVariable = layerData.var<LTSRateAndState::StateVariable>(Cfg(), place);
    data->f0 = layerData.var<LTSRateAndState::RsF0>(Cfg(), place);
    data->muW = layerData.var<LTSRateAndState::RsMuW>(Cfg(), place);
    data->b = layerData.var<LTSRateAndState::RsB>(Cfg(), place);
    data->convergenceInner = layerData.var<LTSRateAndState::ConvergenceInner>(Cfg(), place);
    data->convergenceOuter = layerData.var<LTSRateAndState::ConvergenceOuter>(Cfg(), place);
    LawT::copyLawStorageToLocal(data, layerData, place);
    TPMethod::copyStorageToLocal(data, layerData, place);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void preHook(FrictionLawContext<Cfg, Type>& ctx) {
    // copy state variable from last time step
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      ctx.stateVariableBuffer[pointIndex] = ctx.data->stateVariable[ctx.ltsFace][pointIndex];
    });
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void postHook(FrictionLawContext<Cfg, Type>& ctx) {
    LawT::resampleStateVar(ctx);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void updateFrictionAndSlip(FrictionLawContext<Cfg, Type>& ctx,
                                                       std::uint32_t timeIndex) {
    alignas(pointAlignment<Cfg, Type>()) RateAndStateStep<Cfg, Type> step{};
    // the pressurization of this step of this face
    alignas(pointAlignment<Cfg, Type>()) typename TPMethod::template StepCoefficients<Type>
        tpStep{};

    // compute initial slip rate and reference values
    calcInitialVariables(ctx, step);
    // compute slip rates by solving non-linear system of equations
    updateStateVariableIterative(ctx, step, tpStep, timeIndex);

    // compute final thermal pressure and normalStress
    TPMethod::finalizeFluidPressure(ctx, step.normalStress, step.localSlipRate, timeIndex, tpStep);
    updateDirectionAndProjections(ctx, step);
    updateNormalStress(ctx, step);
    // compute final slip rates and traction from average of the iterative solution and initial
    // guess
    calcSlipRateAndTraction(ctx, step, timeIndex);
  }

  //! the friction coefficient, whose falling to the weakened one is the dynamic stress time
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static PointValues<Cfg, Type>
      dynamicStressQuantity(FrictionLawContext<Cfg, Type>& ctx) {
    PointValues<Cfg, Type> mu{};
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      mu[pointIndex] = ctx.data->mu[ctx.ltsFace][pointIndex];
    });
    return mu;
  }

  /**
   * output time when the friction coefficient has fallen to the weakened one after the rupture
   * arrived, interpolated between the time points of the previous and the current evaluation (see
   * common::crossingTime), and not before the rupture time
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void saveDynamicStressOutput(FrictionLawContext<Cfg, Type>& ctx,
                                                         const PointValues<Cfg, Type>& previousMu,
                                                         const common::SubStepTimes<Cfg>& times) {
    const auto ltsFace = ctx.ltsFace;
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      const real ruptureTime = ctx.data->ruptureTime[ltsFace][pointIndex];
      const real muW = ctx.data->muW[ltsFace][pointIndex];
      const real threshold =
          muW + static_cast<real>(0.05) * (ctx.data->f0[ltsFace][pointIndex] - muW);
      const real mu = ctx.data->mu[ltsFace][pointIndex];
      if (ruptureTime > static_cast<real>(0.0) && ruptureTime <= times.point &&
          ctx.data->dynStressTimePending[ltsFace][pointIndex] && mu <= threshold) {
        ctx.data->dynStressTime[ltsFace][pointIndex] =
            std::max(ruptureTime,
                     common::crossingTime<Cfg>(
                         previousMu[pointIndex], mu, threshold, times.previousPoint, times.point));
        ctx.data->dynStressTimePending[ltsFace][pointIndex] = false;
      }
    });
  }

  /*
   * Compute shear stress magnitude, localSlipRate, effective normal stress, reference state
   * variable. Also sets slipRateMagnitude member to reference value.
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void calcInitialVariables(FrictionLawContext<Cfg, Type>& ctx,
                                                      RateAndStateStep<Cfg, Type>& step) {
    const auto ltsFace = ctx.ltsFace;
    const auto& impAndEta = ctx.data->impAndEta[ltsFace];
    const auto& impedanceMatrices = ctx.data->impedanceMatrices[ltsFace];

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      step.stateVarReference[pointIndex] = ctx.stateVariableBuffer[pointIndex];

      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 =
          ctx.initialStress.traction1[pointIndex] + ctx.faultStresses.traction1[pointIndex];
      const real totalTraction2 =
          ctx.initialStress.traction2[pointIndex] + ctx.faultStresses.traction2[pointIndex];
      const real absoluteTraction = misc::magnitude(totalTraction1, totalTraction2);
      step.absoluteShearTraction[pointIndex] = absoluteTraction;

      const auto [unusedEta, invEta] = common::projectEta<Cfg>(
          impAndEta, impedanceMatrices, totalTraction1, totalTraction2, absoluteTraction);
      step.etaInv[pointIndex] = invEta;

      step.etaNormal[pointIndex] = common::projectEtaNormal<Cfg>(
          impAndEta, impedanceMatrices, totalTraction1, totalTraction2, absoluteTraction);

      // initial slip direction: the trial traction. For isotropy this stays exact.
      const real invAbsolute = (absoluteTraction > 0) ? static_cast<real>(1.0) / absoluteTraction
                                                      : static_cast<real>(0.0);
      step.slipDirection1[pointIndex] = totalTraction1 * invAbsolute;
      step.slipDirection2[pointIndex] = totalTraction2 * invAbsolute;

      // The following process is adapted from that described by Kaneko et al. (2008)
      const real slipRateMagnitude =
          std::max(rs::almostZero<real>(),
                   misc::magnitude(ctx.data->slipRate1[ltsFace][pointIndex],
                                   ctx.data->slipRate2[ltsFace][pointIndex]));
      ctx.data->slipRateMagnitude[ltsFace][pointIndex] = slipRateMagnitude;
      step.localSlipRate[pointIndex] = slipRateMagnitude;
    });

    // after the loop: updateNormalStress reads the slip rate magnitude, which is only set above
    updateNormalStress(ctx, step);
  }

  /**
   * Anisotropy: sweeps the slip direction, and with it every direction dependent projection of the
   * impedance. A no-op for isotropy, where the slip is exactly parallel to the trial traction.
   *
   * The exact condition is tau0 = (S I + V eta_ss) n with |n| = 1, see common::updateSlipDirection.
   * The strength belonging to the current slip rate is recovered as S = n^T tau0 - V * eta_proj,
   * so no additional friction law evaluation is needed here. Riding along with the existing outer
   * fixed-point loop makes the sweep essentially free. It skips the points that `frozen` marks,
   * whose fixed point has converged (see updateStateVariableIterative).
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void updateDirectionAndProjections(
      [[maybe_unused]] FrictionLawContext<Cfg, Type>& ctx,
      [[maybe_unused]] RateAndStateStep<Cfg, Type>& step,
      [[maybe_unused]] const PointArray<std::int32_t, Cfg, Type>* frozen = nullptr) {
    if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
      const auto ltsFace = ctx.ltsFace;
      const auto& impAndEta = ctx.data->impAndEta[ltsFace];
      const auto& impedanceMatrices = ctx.data->impedanceMatrices[ltsFace];
      forEachPoint(ctx, [&](std::uint32_t pointIndex) {
        if (frozen != nullptr && (*frozen)[pointIndex] != 0) {
          return;
        }
        const real totalTraction1 =
            ctx.initialStress.traction1[pointIndex] + ctx.faultStresses.traction1[pointIndex];
        const real totalTraction2 =
            ctx.initialStress.traction2[pointIndex] + ctx.faultStresses.traction2[pointIndex];
        const real trialMagnitude = misc::magnitude(totalTraction1, totalTraction2);
        const real slipRate = ctx.data->slipRateMagnitude[ltsFace][pointIndex];

        const auto [etaProj, unusedInv] = common::projectEta<Cfg>(impAndEta,
                                                                  impedanceMatrices,
                                                                  step.slipDirection1[pointIndex],
                                                                  step.slipDirection2[pointIndex],
                                                                  static_cast<real>(1.0));

        const real strength = step.absoluteShearTraction[pointIndex] - slipRate * etaProj;

        const auto [n1, n2] = common::updateSlipDirection<Cfg>(impAndEta,
                                                               impedanceMatrices,
                                                               strength,
                                                               slipRate,
                                                               totalTraction1,
                                                               totalTraction2,
                                                               trialMagnitude);
        step.slipDirection1[pointIndex] = n1;
        step.slipDirection2[pointIndex] = n2;

        step.absoluteShearTraction[pointIndex] = n1 * totalTraction1 + n2 * totalTraction2;

        const auto [etaUnused, invEta] =
            common::projectEta<Cfg>(impAndEta, impedanceMatrices, n1, n2, static_cast<real>(1.0));
        step.etaInv[pointIndex] = invEta;
        step.etaNormal[pointIndex] = common::projectEtaNormal<Cfg>(
            impAndEta, impedanceMatrices, n1, n2, static_cast<real>(1.0));
      });
    }
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      updateStateVariableIterative(FrictionLawContext<Cfg, Type>& ctx,
                                   RateAndStateStep<Cfg, Type>& step,
                                   TPMethod::template StepCoefficients<Type>& tpStep,
                                   std::uint32_t timeIndex) {
    const auto ltsFace = ctx.ltsFace;
    const real timeIncrement = ctx.args->deltaT[timeIndex];
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> testSlipRate{};
    // whether the fixed point of a point has converged. A converged point stays where it is while
    // the other points of its face (host) or warp (device) iterate on, so that its result does not
    // depend on the points it is evaluated with. int not bool: keeps ICX SIMD happy
    PointArray<std::int32_t, Cfg, Type> convergedOuter{};

    // use:
    // - (inner loop) Newton-Raphson to find the fixed point slip rate with a _fixed_ state and
    // stress
    // - (outer loop) fixed-point iteration to find the fixed point slip rate with varying state and
    // stress; using the previous Newton-Raphson step (why not combine? Mainly because: thermal
    // pressurization might happen which could be a bit expensive to differentiate (though it's
    // doable; just only take all code paths that have tauV in them). Maybe the state can be
    // combined in; however that might need some more nonlinear function evaluations per step, and
    // thus could be slower). The fixed-point iteration is "regularized" by averging it with the
    // previous estimate. I.e. we compute x_(n+1) = (x_n + f(x_n)) / 2. Any fixed point we find is a
    // fixed point of f.

    // procedure source: Kaneko 2008; doi:10.1029/2007JB005154 . Section 2.3. But extended for a
    // virtually unlimited number of outer fixed point iterations.

    // The pressurization of this step is affine in the shear heating, so one walk over the
    // wavenumber grid covers every iteration below; what the iterations then do is evaluate two
    // coefficients per point.
    TPMethod::prepareFluidPressure(ctx, timeIndex, tpStep);

    for (std::uint32_t j = 0; j < ctx.data->drParameters.rsNumberStateVariableUpdates; j++) {
      if constexpr (!LawT::FoldsStateVariable) {
        forEachPoint(ctx, [&](std::uint32_t pointIndex) {
          if (convergedOuter[pointIndex] != 0) {
            return;
          }
          // fault strength using friction coefficient and fluid pressure from previous
          // timestep/iteration update state variable using sliprate from the previous time step
          ctx.stateVariableBuffer[pointIndex] = static_cast<real>(LawT::stateVariableAt(
              ctx,
              pointIndex,
              step.stateVarReference[pointIndex],
              timeIncrement,
              static_cast<LawT::StateScalar>(step.localSlipRate[pointIndex])));
        });
      }
      // the pressure and the normal stress of a converged point stay as well: with the thermal
      // pressurization, they would go on following each other
      TPMethod::applyShearHeating(
          ctx, step.normalStress, step.localSlipRate, tpStep, &convergedOuter);

      updateDirectionAndProjections(ctx, step, &convergedOuter);
      updateNormalStress(ctx, step, &convergedOuter);

      // solve for new slip rate
      invertSlipRateIterative(ctx, step, tpStep, timeIncrement, testSlipRate, convergedOuter);

      forEachPoint(ctx, [&](std::uint32_t pointIndex) {
        if (convergedOuter[pointIndex] != 0) {
          return;
        }
        // update local slip rate, now using V=(Vnew+Vold)/2
        // For the next SV update, use the mean slip rate between the initial guess and the one
        // found (Kaneko 2008, step 6)
        step.localSlipRate[pointIndex] =
            static_cast<real>(0.5) *
            (ctx.data->slipRateMagnitude[ltsFace][pointIndex] + testSlipRate[pointIndex]);

        // Relative, like the criterion of the inner solve: the slip rates this loop walks
        // through span twenty decades, and a fixed step in m/s is a different demand at every
        // one of them -- no requirement at all where a point creeps, and more precision than
        // the inner solve resolves where it slips fast. The tolerance is raised to a few ulp,
        // since no nonzero relative step is smaller.
        const auto stateTolerance =
            std::max(static_cast<real>(ctx.data->drParameters.rsStateTolerance),
                     static_cast<real>(4.0) * std::numeric_limits<real>::epsilon());
        convergedOuter[pointIndex] =
            std::abs(testSlipRate[pointIndex] - ctx.data->slipRateMagnitude[ltsFace][pointIndex]) <=
                    stateTolerance * std::abs(testSlipRate[pointIndex])
                ? 1
                : 0;

        // solve again for Vnew
        ctx.data->slipRateMagnitude[ltsFace][pointIndex] = testSlipRate[pointIndex];
      });

      // exit early (on the device per warp, which prevents the divergence of its threads)
      if (allPoints(ctx, convergedOuter)) {
        break;
      }
    }

    if constexpr (Type == common::RangeType::GPU) {
      deviceBarrier(ctx);
    }

    // update (outer) (non-)convergence
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      ctx.data->convergenceOuter[ltsFace][pointIndex] &= (convergedOuter[pointIndex] != 0);
    });
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void calcSlipRateAndTraction(FrictionLawContext<Cfg, Type>& ctx,
                                                         const RateAndStateStep<Cfg, Type>& step,
                                                         std::uint32_t timeIndex) {
    const auto ltsFace = ctx.ltsFace;
    const real deltaT = ctx.args->deltaT[timeIndex];
    const auto& impAndEta = ctx.data->impAndEta[ltsFace];
    const auto& impedanceMatrices = ctx.data->impedanceMatrices[ltsFace];

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      // SV from mean slip rate in tmp
      ctx.stateVariableBuffer[pointIndex] = static_cast<real>(
          LawT::stateVariableAt(ctx,
                                pointIndex,
                                step.stateVarReference[pointIndex],
                                deltaT,
                                static_cast<LawT::StateScalar>(step.localSlipRate[pointIndex])));
    });

    // Where the law can state its state variable as a function of the slip rate, the friction
    // coefficient follows the accepted slip rate directly and the coefficients that would be
    // precomputed from a supplied state have no meaning here.
    alignas(pointAlignment<Cfg, Type>()) const auto details = frozenMuDetails(ctx);

    forEachPoint<rs_detail::VectorizeStrict>(ctx, [&](std::uint32_t pointIndex) {
      const real slipRateMagnitude = ctx.data->slipRateMagnitude[ltsFace][pointIndex];
      // update LocMu for next strength determination, only needed for last update
      real mu{};
      if constexpr (LawT::FoldsStateVariable) {
        mu = LawT::updateMuFolded(
            ctx, pointIndex, slipRateMagnitude, step.stateVarReference[pointIndex], deltaT);
      } else {
        mu = LawT::updateMu(ctx, pointIndex, slipRateMagnitude, details);
      }
      ctx.data->mu[ltsFace][pointIndex] = mu;
      const real strength = -mu * step.normalStress[pointIndex];

      // the direction along which the slip rate is decomposed; scaled such that dividing by
      // `divisor` yields the unit slip direction. For isotropy slipDirection is the (normalized)
      // trial traction and absoluteShearTraction is its magnitude.
      const real dirTraction1 =
          step.slipDirection1[pointIndex] * step.absoluteShearTraction[pointIndex];
      const real dirTraction2 =
          step.slipDirection2[pointIndex] * step.absoluteShearTraction[pointIndex];

      const auto [eta, unusedInv] = common::projectEta<Cfg>(impAndEta,
                                                            impedanceMatrices,
                                                            step.slipDirection1[pointIndex],
                                                            step.slipDirection2[pointIndex],
                                                            static_cast<real>(1.0));

      const auto divisor = strength + eta * slipRateMagnitude;
      const real slipRate1 = slipRateMagnitude * dirTraction1 / divisor;
      const real slipRate2 = slipRateMagnitude * dirTraction2 / divisor;
      ctx.data->slipRate1[ltsFace][pointIndex] = slipRate1;
      ctx.data->slipRate2[ltsFace][pointIndex] = slipRate2;

      const auto [tU1, tU2] =
          common::matmulEta<Cfg>(impAndEta, impedanceMatrices, slipRate1, slipRate2);

      // calculate traction
      // note that the normal stress written here is the *dynamic* normal traction, i.e. in the
      // same space as faultStresses/qInterpolated -- not the effective normal stress used for the
      // friction strength above (which additionally carries the initial stress, the initial
      // pressure, thermal pressurization and the min(., 0) clamp).
      ctx.tractionResults.normalStress[pointIndex] = ctx.faultStresses.normalStress[pointIndex] -
                                                     slipRateMagnitude * step.etaNormal[pointIndex];
      ctx.tractionResults.traction1[pointIndex] = ctx.faultStresses.traction1[pointIndex] - tU1;
      ctx.tractionResults.traction2[pointIndex] = ctx.faultStresses.traction2[pointIndex] - tU2;
      ctx.data->traction1[ltsFace][pointIndex] = ctx.tractionResults.traction1[pointIndex];
      ctx.data->traction2[ltsFace][pointIndex] = ctx.tractionResults.traction2[pointIndex];

      // Compute slip
      // ABS of locSlipRate removed as it would be the accumulated slip that is usually not needed
      // in the solver, see linear slip weakening
      ctx.data->accumulatedSlipMagnitude[ltsFace][pointIndex] += slipRateMagnitude * deltaT;

      // update directional slip
      ctx.data->slip1[ltsFace][pointIndex] += slipRate1 * deltaT;
      ctx.data->slip2[ltsFace][pointIndex] += slipRate2 * deltaT;
    });
  }

  /**
   * Solve for new slip rate (\f$\hat{s}\f$) with a bracketed Newton (rtsafe). \f$\hat{s}\f$ solves
   * \f[g := -\frac{1}{\eta}(|\sigma(\hat{s})|\,\mu(\hat{s}) - \Theta) - \hat{s} = 0,\f]
   * c.f. Carsten Uphoff's dissertation eq. (4.57), with \f$\eta\f$ projected onto the fault plane
   * per point and \f$\sigma\f$ following the slip rate through the anisotropic normal coupling.
   * The root is bracketed in closed form, and the bracket survives that coupling because it only
   * needs \f$\mu(0)=0\f$, \f$\mu\ge0\f$ and \f$|\sigma|\ge0\f$ (endpoints are NOT evaluated):
   *   g(0+)           = +invEta * Theta            > 0
   *   g(Theta*invEta) = -invEta * |sigma| * mu    <= 0.
   * The root need not be unique, and for a law that folds its state variable it usually is not. The
   * state variable weakens with the slip rate, so once it travels inside the residual
   * \f$\partial\mu/\partial V\f$ changes sign -- around a metre per second for the fast law --
   * and \f$g'\f$ rises towards zero with it: a plain isotropic fault at 50 MPa reaches
   * \f$g' = -0.72\f$ there, and stepping rupture histories through the fast law leaves 471 of 3240
   * inversions with more than one root inside the bracket. The anisotropic normal coupling weakens
   * \f$g'\f$ further (cf. the derivative below). Bisection converges to a root in the bracket
   * either way, so the solver does not rest on uniqueness -- but which root it converges to is
   * decided by where it starts, which is why the warm start below matters for more than speed:
   * started at the previous step 2 of those 3240 land on a different root than a bisection from the
   * whole bracket would, started at the free-slip limit 469 of them do.
   * We take Newton while it stays in the bracket and outruns bisection, else bisect. The bracket
   * is non-increasing and loses half of its decades on every fallback, so the iterate settles and
   * termination is relative in SLIP-RATE space (|dV| < xacc * V). Two floors keep that test
   * reachable in finite precision: xacc is clamped to a few ulp, and an iterate whose residual has
   * sunk into the rounding noise of its own evaluation counts as converged, because no further step
   * can be told from noise.
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      invertSlipRateIterative(FrictionLawContext<Cfg, Type>& ctx,
                              const RateAndStateStep<Cfg, Type>& step,
                              const TPMethod::template StepCoefficients<Type>& tpStep,
                              real timeIncrement,
                              PointValues<Cfg, Type>& slipRateTest,
                              const PointArray<std::int32_t, Cfg, Type>& frozen) {
    const auto ltsFace = ctx.ltsFace;

    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> muF{};
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> g{};
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> xLow{};
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> xHigh{};
    // previous step, for the "outrun bisection" test
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> dxOld{};
    // rounding noise of the residual, per point
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> gNoise{};
    // d(mu)/dV, from the same pass as mu itself
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> dMuF{};
    // |sigma| at the trial slip rate
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> absSigma{};
    // and its derivative, from the same pass
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> dAbsSigma{};
    // int not bool: keeps ICX SIMD happy (cf. below)
    PointArray<std::int32_t, Cfg, Type> converged{};

    // Number of roundings that enter one residual evaluation; used to size both floors below.
    constexpr real NoiseFactor = 4;
    constexpr real Eps = std::numeric_limits<real>::epsilon();

    // rsSlipRateTolerance is a RELATIVE step tolerance. A nonzero step is at least one ulp of the
    // iterate, so a tolerance below Eps cannot be met at all and would leave the exact-fixed-point
    // guard as the only way out. Clamp it to a few ulp.
    const real xacc =
        std::max(static_cast<real>(ctx.data->drParameters.rsSlipRateTolerance), NoiseFactor * Eps);

    // Where the law can state its state variable as a function of the slip rate, the residual
    // evaluates it at the trial slip rate rather than at one the outer fixed point supplies, and
    // the coefficients that would be precomputed from a frozen state have no meaning here.
    alignas(pointAlignment<Cfg, Type>()) const auto details = frozenMuDetails(ctx);

    // closed-form bracket + warm start (clamped previous-step V); no endpoint evaluations
    forEachPoint<rs_detail::VectorizeStrict>(ctx, [&](std::uint32_t pointIndex) {
      const real lo = rs::almostZero<real>();
      const real hi =
          std::max(lo, step.absoluteShearTraction[pointIndex] * step.etaInv[pointIndex]);
      // A point that carries no normal stress at the free-slip limit has its root exactly there:
      // |sigma| vanishes, g is the line Theta * invEta - V, and g(hi) = 0. That is worth taking
      // directly, because it is the one root rtsafe cannot approach: a root on the bracket
      // boundary leaves the Newton step the same size as the previous one, so the guard falls
      // back to bisection on every iteration and the solve spends its whole budget halving.
      const bool openAtLimit =
          stickAt(step, TPMethod::fluidPressureOffset(tpStep, pointIndex), hi, pointIndex) >=
          static_cast<real>(0.0);
      xLow[pointIndex] = lo;
      xHigh[pointIndex] = hi;
      slipRateTest[pointIndex] =
          openAtLimit
              ? hi
              : std::min(std::max(ctx.data->slipRateMagnitude[ltsFace][pointIndex], lo), hi);
      dxOld[pointIndex] = hi - lo;
      // a point whose outer fixed point has converged is not solved again
      converged[pointIndex] = (openAtLimit || frozen[pointIndex] != 0) ? 1 : 0;
    });

    for (std::uint32_t i = 0; i < ctx.data->drParameters.rsMaxNumberSlipRateUpdates; i++) {

      // residual + bracket maintenance. Transcendental mu() eval stays vectorized and unmasked;
      // converged lanes' brackets are simply never read again, so no mask is needed here.
      // >>> precision knob: evaluate muF/g (and dG below) in double -- promote |sigma| and Theta
      //     (lossless float->double) -- to drop the noise floor AND make the sign(g) test exact.
      forEachPoint<rs_detail::VectorizeStrict>(ctx, [&](std::uint32_t pointIndex) {
        const real x = slipRateTest[pointIndex];
        // one pass through mu() yields the value and its derivative; the friction law is written
        // once and instantiated for a dual number here. Where the state variable is folded in, the
        // derivative that comes back is the one of the composition mu(V, psi(V)).
        const Dual<real> trial(x, static_cast<real>(1.0));
        Dual<real> mu{};
        if constexpr (LawT::FoldsStateVariable) {
          mu = LawT::updateMuFolded(
              ctx, pointIndex, trial, step.stateVarReference[pointIndex], timeIncrement);
        } else {
          mu = LawT::updateMu(ctx, pointIndex, trial, details);
        }
        muF[pointIndex] = mu.value;
        dMuF[pointIndex] = mu.derivative;
        // sigma follows the trial slip rate, so it is evaluated at x rather than taken frozen:
        // that moves the normal coupling and the pressurization out of the outer fixed point and
        // into this Newton. The dual carries d|sigma|/dV of the whole composition, including the
        // part that reaches through mu.
        const Dual<real> sigma = rs::effectiveNormalStress(
            stickAt(step, TPMethod::fluidPressureOffset(tpStep, pointIndex), trial, pointIndex),
            TPMethod::fluidPressureSlope(tpStep, pointIndex),
            trial,
            mu);
        const auto absSigmaDual = abs(sigma);
        absSigma[pointIndex] = absSigmaDual.value;
        dAbsSigma[pointIndex] = absSigmaDual.derivative;
        g[pointIndex] = -step.etaInv[pointIndex] * (absSigma[pointIndex] * muF[pointIndex] -
                                                    step.absoluteShearTraction[pointIndex]) -
                        x;
        // |sigma| * mu and tau cancel at the root, so the rounding error of g does not shrink
        // with the iterate: it stays at Eps times the magnitude of the two cancelling terms. Below
        // that level the sign of g -- and with it the bracket update -- carries no information.
        gNoise[pointIndex] =
            NoiseFactor * Eps * step.etaInv[pointIndex] *
            (absSigma[pointIndex] * muF[pointIndex] + step.absoluteShearTraction[pointIndex]);
        // g decreasing: g>0 => root above x
        const bool gPos = g[pointIndex] > static_cast<real>(0);
        xLow[pointIndex] = gPos ? x : xLow[pointIndex];
        xHigh[pointIndex] = !gPos ? x : xHigh[pointIndex];
      });

      // Newton-or-bisect select (contains the division: guarded by the *non-strict* ICX macro).
      forEachPoint<rs_detail::VectorizeSelect>(ctx, [&](std::uint32_t pointIndex) {
        const real x = slipRateTest[pointIndex];

        const real dGFrozen = -step.etaInv[pointIndex] * (absSigma[pointIndex] * dMuF[pointIndex]) -
                              static_cast<real>(1.0);
        const real dGCoupled =
            -step.etaInv[pointIndex] * (absSigma[pointIndex] * dMuF[pointIndex] +
                                        dAbsSigma[pointIndex] * muF[pointIndex]) -
            static_cast<real>(1.0);
        // A fault that loses normal stress as it slips (etaNormal < 0) is the only case in which
        // the coupling can weaken g. It stays strictly decreasing as long as
        //   |etaNormal| * mu < eta_proj + |sigma| * mu' ,
        // which is NOT a comfortable margin: scanning dGCoupled over the whole bracket, a locked
        // point at |etaNormal| = eta_proj already peaks at dGCoupled = -0.06, and turns positive
        // slightly above that ratio; a slipping or nearly open point holds out to roughly three
        // times eta_proj. Positive definiteness of the 3x3 impedance bounds the ratio by
        // sqrt(eta_nn / eta_proj), so a strongly anisotropic material lands close to the edge.
        // Past it g has several roots in the bracket and the solve returns one of them -- still a
        // solution of the discrete problem, but not necessarily the intended branch. This is why
        // the bracket, rather than the derivative, carries the robustness here.
        // For etaNormal > 0 the coupled derivative is the better conditioned of the two. dGFrozen
        // is negative by construction, which keeps the Newton step pointing into the bracket even
        // when the coupled one would not.
        const real dG = (dGCoupled < static_cast<real>(0.0)) ? dGCoupled : dGFrozen;

        // Which variable the residual is closer to linear in is written in dG itself: it carries
        // the strength term and the -1 of the radiation damping. Where the strength term wins,
        // mu is in its logarithmic branch and g is nearly linear in log V; where the damping
        // wins, g is nearly linear in V. The two are equal at dG = -2.
        //
        // The log step is Newton on h(u) = g(exp(u)), h'(u) = dG * V, hence multiplicative. The
        // clamp keeps exp() inside the range before the bracket test gets to reject the step, and
        // a multiplicative step cannot leave the positive axis.
        real xNewton{};
        if (dG < static_cast<real>(-2.0)) {
          const real du = std::min(std::max(-g[pointIndex] / (dG * x), static_cast<real>(-60.0)),
                                   static_cast<real>(60.0));
          xNewton = x * std::exp(du);
        } else {
          xNewton = x - g[pointIndex] / dG;
        }
        // Bisect geometrically. The bracket spans the whole admissible range of slip rates,
        // from almostZero() up to the free-slip limit tau/eta_s, so its arithmetic midpoint sits
        // many orders of magnitude above the root of a locked or creeping point, and a fallback
        // would then need one halving per factor of two to walk back down. The geometric midpoint
        // halves the number of decades instead, which is the scale the root lives on. The two
        // square roots keep the product from underflowing for the smallest brackets.
        const real xBisect = std::sqrt(xLow[pointIndex]) * std::sqrt(xHigh[pointIndex]);
        const bool useBisect =
            (xNewton <= xLow[pointIndex]) || (xNewton >= xHigh[pointIndex]) ||
            (std::abs(static_cast<real>(2.0) * g[pointIndex]) > std::abs(dxOld[pointIndex] * dG));
        const real xUpdated = useBisect ? xBisect : xNewton;
        const real stepSize = xUpdated - x;

        // V-space convergence, residual-noise floor, and the no-representable-change guard
        // => cannot livelock
        const bool laneConv = (std::abs(stepSize) < xacc * std::abs(x)) ||
                              (std::abs(g[pointIndex]) <= gNoise[pointIndex]) || (xUpdated == x);
        const std::int32_t nowConv = (converged[pointIndex] != 0 || laneConv) ? 1 : 0;
        converged[pointIndex] = nowConv;

        // freeze converged lanes at the accepted point (so muF stays the consistent pair);
        // advance the rest
        slipRateTest[pointIndex] = (nowConv != 0) ? x : xUpdated;
        dxOld[pointIndex] = stepSize;
      });

      // exit early (on the device per warp, which prevents the divergence of its threads)
      if (allPoints(ctx, converged)) {
        break;
      }
    }

    // publish mu consistent with the accepted slip rate
    // in case of non-convergence, flag the offending lanes -- but now
    // keyed on the x-criterion, consistent with the termination test above
    forEachPoint<rs_detail::VectorizeStrict>(ctx, [&](std::uint32_t pointIndex) {
      if (frozen[pointIndex] != 0) {
        return;
      }
      if constexpr (LawT::FoldsStateVariable) {
        ctx.data->mu[ltsFace][pointIndex] = LawT::updateMuFolded(ctx,
                                                                 pointIndex,
                                                                 slipRateTest[pointIndex],
                                                                 step.stateVarReference[pointIndex],
                                                                 timeIncrement);
      } else {
        ctx.data->mu[ltsFace][pointIndex] =
            LawT::updateMu(ctx, pointIndex, slipRateTest[pointIndex], details);
      }
      ctx.data->convergenceInner[ltsFace][pointIndex] &= (converged[pointIndex] != 0);
    });
  }

  /// the coefficients of mu precomputed from the state variable of the points; empty for a law that
  /// folds its state variable into the residual, as there is no frozen one to precompute them from
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static auto frozenMuDetails(const FrictionLawContext<Cfg, Type>& ctx) {
    if constexpr (!LawT::FoldsStateVariable) {
      return LawT::getMuDetails(ctx);
    } else {
      return typename LawT::template MuDetails<Type>{};
    }
  }

  /**
   * Effective normal stress, including the anisotropic normal/shear coupling.
   *
   * With the 3x3 impedance eta, shear slip changes the fault-normal traction:
   *   sigma(V) = sigma_stick - V * etaNormal,   etaNormal = (eta * n)_n
   * For every isotropic material etaNormal is zero and this reduces to the previous formula.
   *
   * The slip rate is taken from slipRateMagnitude, i.e. from the previous outer fixed-point
   * iteration (or, on entry, from the previous time step). normalStressStick keeps the part that
   * does not depend on it, so that the Newton solve can follow sigma(V) itself. The points that
   * `frozen` marks keep theirs (see updateStateVariableIterative).
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      updateNormalStress(FrictionLawContext<Cfg, Type>& ctx,
                         RateAndStateStep<Cfg, Type>& step,
                         const PointArray<std::int32_t, Cfg, Type>* frozen = nullptr) {
    const auto ltsFace = ctx.ltsFace;
    // Todo(SW): consider poroelastic materials together with thermal pressurization
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      if (frozen != nullptr && (*frozen)[pointIndex] != 0) {
        return;
      }
      // the stick keeps what does not follow this step's shear heating; the pressurization
      // enters where sigma is formed, so that the solve can resolve it instead of lagging it
      step.normalStressStick[pointIndex] =
          ctx.faultStresses.normalStress[pointIndex] + ctx.initialStress.normalStress[pointIndex] +
          ctx.faultStresses.fluidPressure[pointIndex] + ctx.initialStress.fluidPressure[pointIndex];
      step.normalStress[pointIndex] = std::min(
          static_cast<real>(0.0),
          step.normalStressStick[pointIndex] - TPMethod::getFluidPressure(ctx, pointIndex) -
              ctx.data->slipRateMagnitude[ltsFace][pointIndex] * step.etaNormal[pointIndex]);
    });
  }

  /**
   * The normal stress a trial slip rate belongs to, before the pressurization of this step: the
   * stick, less the pore pressure the step reaches without shear heating (pressureOffset), and less
   * the normal traction the slip moves through the anisotropic coupling, which vanishes for every
   * other material. rs::effectiveNormalStress adds the pressurization by the shear heating of the
   * trial slip rate. The effective normal stress of updateNormalStress is no stick here: it already
   * holds the pressure of the shear heating of the previous iterate.
   */
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S stickAt(const RateAndStateStep<Cfg, Type>& step,
                                      real pressureOffset,
                                      S slipRate,
                                      std::uint32_t pointIndex) {
    return S(step.normalStressStick[pointIndex] - pressureOffset) -
           slipRate * S(step.etaNormal[pointIndex]);
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_RATEANDSTATE_H_
