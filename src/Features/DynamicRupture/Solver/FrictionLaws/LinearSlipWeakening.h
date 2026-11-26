// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_LINEARSLIPWEAKENING_H_

#include "Core/Common/Marker.h"
#include "Core/Common/MaterialType.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

/**
 * The linear slip weakening friction laws, once for the host and the device: the hooks of a law for
 * FrictionLawDriver, which the solvers of the host (cpu::LinearSlipWeakeningLaw) and of the device
 * (gpu::LinearSlipWeakeningLaw) run. The specialization `SpecializationT` adapts the strength, the
 * state variable and the resampling of the slip rate (NoSpecialization, BiMaterialFault, TPApprox).
 */
namespace seissol::dr::friction_law {

template <typename Cfg, typename SpecializationT>
class LinearSlipWeakening {
  public:
  using real = Real<Cfg>;

  // the slip rate of a sub-step is the one the friction solve finds at its time point
  static constexpr bool SlipRateAtIntervalEnds = false;

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData,
                                             seissol::initializer::AllocationPlace place) {
    data->dC = layerData.var<LTSLinearSlipWeakening::DC>(Cfg(), place);
    data->muS = layerData.var<LTSLinearSlipWeakening::MuS>(Cfg(), place);
    data->muD = layerData.var<LTSLinearSlipWeakening::MuD>(Cfg(), place);
    data->cohesion = layerData.var<LTSLinearSlipWeakening::Cohesion>(Cfg(), place);
    data->forcedRuptureTime =
        layerData.var<LTSLinearSlipWeakening::ForcedRuptureTime>(Cfg(), place);
    SpecializationT::copyStorageToLocal(data, layerData, place);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void preHook(FrictionLawContext<Cfg, Type>& /*ctx*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void postHook(FrictionLawContext<Cfg, Type>& /*ctx*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void updateFrictionAndSlip(FrictionLawContext<Cfg, Type>& ctx,
                                                       uint32_t timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    calcStrengthHook(ctx, timeIndex);
    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    calcSlipRateAndTraction(ctx, timeIndex);
    // integrate state variable in time
    calcStateVariableHook(ctx, timeIndex);
    // compute friction coefficient based on state variable and slip
    frictionFunctionHook(ctx);
  }

  /**
   * calculate fault strength (Uphoff eq 2.44) with addition cohesion term.
   * The anisotropic normal/shear coupling is deliberately not applied here: since the strength is
   * affine in the normal stress, it is handled exactly through the divisor in
   * calcSlipRateAndTraction, using the slope filled in below.
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void calcStrengthHook(FrictionLawContext<Cfg, Type>& ctx,
                                                  uint32_t timeIndex) {
    const auto deltaT{ctx.args->deltaT[timeIndex]};
    const auto vStar{ctx.data->drParameters.vStar};
    const auto prakashLength{ctx.data->drParameters.prakashLength};
    const auto ltsFace = ctx.ltsFace;

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      auto& strength = ctx.strengthBuffer[pointIndex];

      const real totalNormalStress =
          ctx.initialStress.normalStress[pointIndex] + ctx.faultStresses.normalStress[pointIndex] +
          ctx.initialStress.fluidPressure[pointIndex] + ctx.faultStresses.fluidPressure[pointIndex];

      strength =
          -ctx.data->cohesion[ltsFace][pointIndex] -
          ctx.data->mu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<real>(0.0));

      strength = SpecializationT::strengthHook(ctx,
                                               strength,
                                               ctx.data->slipRateMagnitude[ltsFace][pointIndex],
                                               deltaT,
                                               vStar,
                                               prakashLength,
                                               pointIndex);

      if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
        // d(strength) / d(-sigma_eff). Zero while the normal stress is clamped -- note that the
        // clamp is evaluated at the uncorrected normal stress, which is second order in the
        // coupling and only matters right at the opening threshold.
        ctx.strengthSlopeBuffer[pointIndex] =
            (totalNormalStress < 0 ? ctx.data->mu[ltsFace][pointIndex] : static_cast<real>(0.0)) *
            SpecializationT::strengthHookSlope(
                ctx.data->slipRateMagnitude[ltsFace][pointIndex], deltaT, vStar, prakashLength);
      }
    });
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void calcSlipRateAndTraction(FrictionLawContext<Cfg, Type>& ctx,
                                                         uint32_t timeIndex) {
    const auto deltaT{ctx.args->deltaT[timeIndex]};
    const auto ltsFace = ctx.ltsFace;
    const auto& impAndEta = ctx.data->impAndEta[ltsFace];
    const auto& impedanceMatrices = ctx.data->impedanceMatrices[ltsFace];

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      const real strength = ctx.strengthBuffer[pointIndex];

      // calculate absolute value of stress in Y and Z direction
      const real totalTraction1 =
          ctx.initialStress.traction1[pointIndex] + ctx.faultStresses.traction1[pointIndex];
      const real totalTraction2 =
          ctx.initialStress.traction2[pointIndex] + ctx.faultStresses.traction2[pointIndex];
      const real absoluteTraction = misc::magnitude(totalTraction1, totalTraction2);

      const auto [eta, invEta] = common::projectEta<Cfg>(
          impAndEta, impedanceMatrices, totalTraction1, totalTraction2, absoluteTraction);

      // the direction along which the slip rate is decomposed further down; scaled such that
      // dividing by `divisor` yields the unit slip direction. Only the anisotropic branch below
      // changes them, which the check does not see in the other materials.
      // NOLINTBEGIN(misc-const-correctness)
      real dirTraction1 = totalTraction1;
      real dirTraction2 = totalTraction2;
      real etaEff = eta;
      // NOLINTEND(misc-const-correctness)
      real slipRateMagnitude{};

      if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
        common::SlipRateSolution<Cfg> solution{};
        common::solveSlipRate<Cfg>(impAndEta,
                                   impedanceMatrices,
                                   totalTraction1,
                                   totalTraction2,
                                   absoluteTraction,
                                   strength,
                                   ctx.strengthSlopeBuffer[pointIndex],
                                   solution);
        slipRateMagnitude = solution.slipRate;
        etaEff = solution.etaEff;

        // divisor below equals projectedTraction, so this restores V * n
        dirTraction1 = solution.direction1 * solution.projectedTraction;
        dirTraction2 = solution.direction2 * solution.projectedTraction;
      } else {
        slipRateMagnitude =
            std::max(static_cast<real>(0.0), (absoluteTraction - strength) * invEta);
      }

      // calculate slip rates
      ctx.data->slipRateMagnitude[ltsFace][pointIndex] = slipRateMagnitude;

      const auto divisor = strength + etaEff * slipRateMagnitude;
      ctx.data->slipRate1[ltsFace][pointIndex] = slipRateMagnitude * dirTraction1 / divisor;
      ctx.data->slipRate2[ltsFace][pointIndex] = slipRateMagnitude * dirTraction2 / divisor;

      const auto [tU1, tU2] = common::matmulEta<Cfg>(impAndEta,
                                                     impedanceMatrices,
                                                     ctx.data->slipRate1[ltsFace][pointIndex],
                                                     ctx.data->slipRate2[ltsFace][pointIndex]);

      const auto tUN = common::matmulEtaNormal<Cfg>(impAndEta,
                                                    impedanceMatrices,
                                                    ctx.data->slipRate1[ltsFace][pointIndex],
                                                    ctx.data->slipRate2[ltsFace][pointIndex]);

      // calculate traction
      // the normal stress written here is the *dynamic* normal traction, i.e. in the same space as
      // faultStresses/qInterpolated -- unlike the effective normal stress used for the strength,
      // it carries neither the initial stress nor the min(., 0) clamp.
      ctx.tractionResults.normalStress[pointIndex] =
          ctx.faultStresses.normalStress[pointIndex] - tUN;
      ctx.tractionResults.traction1[pointIndex] = ctx.faultStresses.traction1[pointIndex] - tU1;
      ctx.tractionResults.traction2[pointIndex] = ctx.faultStresses.traction2[pointIndex] - tU2;
      ctx.data->traction1[ltsFace][pointIndex] = ctx.tractionResults.traction1[pointIndex];
      ctx.data->traction2[ltsFace][pointIndex] = ctx.tractionResults.traction2[pointIndex];

      // update directional slip
      ctx.data->slip1[ltsFace][pointIndex] += ctx.data->slipRate1[ltsFace][pointIndex] * deltaT;
      ctx.data->slip2[ltsFace][pointIndex] += ctx.data->slipRate2[ltsFace][pointIndex] * deltaT;
    });
  }

  /**
   * integrates the slip rate in time to the slip, which is the state variable of the law (divided
   * by the critical distance), and lets the forced rupture time weaken the fault as well
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void calcStateVariableHook(FrictionLawContext<Cfg, Type>& ctx,
                                                       uint32_t timeIndex) {
    const auto deltaT{ctx.args->deltaT[timeIndex]};
    const auto riseTime{ctx.data->drParameters.forcedRuptureRiseTime};
    const auto tpProxyExponent{ctx.data->drParameters.tpProxyExponent};
    const auto ltsFace = ctx.ltsFace;

    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> resampledSlipRate{};
    SpecializationT::resampleSlipRate(ctx, ctx.data->slipRateMagnitude[ltsFace], resampledSlipRate);

    real time = ctx.fullUpdateTime;
    for (uint32_t i = 0; i <= timeIndex; ++i) {
      time += ctx.args->deltaT[i];
    }

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      auto& stateVariable = ctx.stateVariableBuffer[pointIndex];

      // integrate slip rate to get slip = state variable
      const auto update = resampledSlipRate[pointIndex] * deltaT;
      ctx.data->accumulatedSlipMagnitude[ltsFace][pointIndex] += update;

      // Actually slip is already the stateVariable for this FL, but to simplify the next equations
      // we divide it here by the critical distance.
      stateVariable = SpecializationT::stateVariableHook(
          ctx.data->accumulatedSlipMagnitude[ltsFace][pointIndex],
          ctx.data->dC[ltsFace][pointIndex],
          tpProxyExponent);

      // Forced rupture time
      real f2 = 0.0;
      if (riseTime == 0) {
        // avoid branching
        // if time > forcedRuptureTime, then f2 = 1.0, else f2 = 0.0
        f2 = static_cast<real>(time >= ctx.data->forcedRuptureTime[ltsFace][pointIndex]);
      } else {
        f2 = misc::clamp((time - ctx.data->forcedRuptureTime[ltsFace][pointIndex]) / riseTime,
                         static_cast<real>(0.0),
                         static_cast<real>(1.0));
      }
      stateVariable = std::max(stateVariable, f2);
    });
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void frictionFunctionHook(FrictionLawContext<Cfg, Type>& ctx) {
    const auto healingThreshold{ctx.data->drParameters.healingThreshold};
    const auto ltsFace = ctx.ltsFace;

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      auto& stateVariable = ctx.stateVariableBuffer[pointIndex];
      ctx.data->mu[ltsFace][pointIndex] =
          ctx.data->muS[ltsFace][pointIndex] -
          (ctx.data->muS[ltsFace][pointIndex] - ctx.data->muD[ltsFace][pointIndex]) * stateVariable;
      // instantaneous healing
      if ((ctx.data->peakSlipRate[ltsFace][pointIndex] > healingThreshold) &&
          (ctx.data->slipRateMagnitude[ltsFace][pointIndex] < healingThreshold)) {
        ctx.data->mu[ltsFace][pointIndex] = ctx.data->muS[ltsFace][pointIndex];
        stateVariable = 0.0;
      }
    });
  }

  //! the slip, whose reaching d_c is the dynamic stress time
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static PointValues<Cfg, Type>
      dynamicStressQuantity(FrictionLawContext<Cfg, Type>& ctx) {
    PointValues<Cfg, Type> slip{};
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      slip[pointIndex] = ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][pointIndex];
    });
    return slip;
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived: when the
   * slip reaches d_c, interpolated over the interval the sub-step integrates the slip over (see
   * common::crossingTime)
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void saveDynamicStressOutput(FrictionLawContext<Cfg, Type>& ctx,
                                                         const PointValues<Cfg, Type>& previousSlip,
                                                         const common::SubStepTimes<Cfg>& times) {
    const auto ltsFace = ctx.ltsFace;

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      const real slip = std::fabs(ctx.data->accumulatedSlipMagnitude[ltsFace][pointIndex]);
      const real dC = ctx.data->dC[ltsFace][pointIndex];
      if (ctx.data->dynStressTimePending[ltsFace][pointIndex] && slip >= dC) {
        ctx.data->dynStressTime[ltsFace][pointIndex] = common::crossingTime<Cfg>(
            std::fabs(previousSlip[pointIndex]), slip, dC, times.start, times.end);
        ctx.data->dynStressTimePending[ltsFace][pointIndex] = false;
      }
    });
  }
};

/**
 * The linear slip weakening law itself: the slip rate is resampled into the polynomial space of
 * the unknowns, so that the slip it integrates to lies there as well.
 */
template <typename Cfg>
class NoSpecialization {
  public:
  using real = Real<Cfg>;

  static void copyStorageToLocal(FrictionLawData<Cfg>* /*data*/,
                                 DynamicRupture::Layer& /*layerData*/,
                                 seissol::initializer::AllocationPlace /*place*/) {}

  /**
   * Resample slip-rate, such that the state increment (slip) lies in the same polynomial space as
   * the degrees of freedom resampleMatrix first projects LocSR on the two-dimensional basis on
   * the reference triangle with degree less or equal than ConvergenceOrder-1, and then evaluates
   * the polynomial at the quadrature points
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      resampleSlipRate(FrictionLawContext<Cfg, Type>& ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints<Cfg>],
                       PointValues<Cfg, Type>& resampledSlipRate) {
    resample(ctx, slipRateMagnitude, resampledSlipRate);
  }

  SEISSOL_HOSTDEVICE static real
      stateVariableHook(real localAccumulatedSlip, real localDc, real /*tpProxyExponent*/) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real strengthHook(FrictionLawContext<Cfg, Type>& /*ctx*/,
                                              real strength,
                                              real /*localSlipRate*/,
                                              real /*deltaT*/,
                                              real /*vStar*/,
                                              real /*prakashLength*/,
                                              std::uint32_t /*pointIndex*/) {
    return strength;
  }

  /**
   * d(strengthHook output) / d(its faultStrength argument).
   *
   * Only needed for the anisotropic normal/shear coupling, where the fault strength depends on the
   * slip rate through the normal stress. MUST be free of side effects -- unlike strengthHook,
   * which advances the Prakash-Clifton state -- and has to be evaluated with the same
   * localSlipRate and deltaT as the corresponding strengthHook call.
   */
  SEISSOL_HOSTDEVICE static real strengthHookSlope(real /*localSlipRate*/,
                                                   real /*deltaT*/,
                                                   real /*vStar*/,
                                                   real /*prakashLength*/) {
    return static_cast<real>(1.0);
  }
};

/**
 * Law for bimaterial faults, implements strength regularization (according to Prakash-Clifton)
 */
template <typename Cfg>
class BiMaterialFault {
  public:
  using real = Real<Cfg>;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data,
                                 DynamicRupture::Layer& layerData,
                                 seissol::initializer::AllocationPlace place) {
    data->regularizedStrength =
        layerData.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(Cfg(), place);
  }

  /*
   * Resampling of the sliprate introduces artificial oscillations into the solution, if we use it
   * together with Prakash-Clifton regularization, so for the BiMaterialFault specialization, we
   * replace the resampling with a simple copy.
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      resampleSlipRate(FrictionLawContext<Cfg, Type>& ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints<Cfg>],
                       PointValues<Cfg, Type>& resampledSlipRate) {
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      resampledSlipRate[pointIndex] = slipRateMagnitude[pointIndex];
    });
  }

  SEISSOL_HOSTDEVICE static real
      stateVariableHook(real localAccumulatedSlip, real localDc, real /*tpProxyExponent*/) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  }

  /// modify strength according to Prakash-Clifton
  /// see e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real strengthHook(FrictionLawContext<Cfg, Type>& ctx,
                                              real faultStrength,
                                              real localSlipRate,
                                              real deltaT,
                                              real vStar,
                                              real prakashLength,
                                              std::uint32_t pointIndex) {
    const auto expval =
        -(std::max(static_cast<real>(0.0), localSlipRate) + vStar) * deltaT / prakashLength;
    const real exp1mterm = -std::expm1(expval);

    // the relaxation towards faultStrength with expm1 alone, which keeps faultStrength as its exact
    // fixed point in any precision (see FastVelocityWeakeningLaw::updateStateVariable)
    const real regularizedStrength = ctx.data->regularizedStrength[ctx.ltsFace][pointIndex];
    const real newStrength =
        regularizedStrength + (faultStrength - regularizedStrength) * exp1mterm;

    ctx.data->regularizedStrength[ctx.ltsFace][pointIndex] = newStrength;
    return newStrength;
  }

  /**
   * See NoSpecialization::strengthHookSlope. The Prakash-Clifton regularization low-passes the
   * strength, so only the fraction exp1mterm of a change in faultStrength arrives instantaneously.
   * regularizedStrength carries the previous step and does not depend on the current normal
   * stress, hence it drops out of the derivative.
   */
  SEISSOL_HOSTDEVICE static real
      strengthHookSlope(real localSlipRate, real deltaT, real vStar, real prakashLength) {
    const auto expval =
        -(std::max(static_cast<real>(0.0), localSlipRate) + vStar) * deltaT / prakashLength;
    return -std::expm1(expval);
  }
};

/**
 * Modified LSW friction as discussed in github issue #1058
 */
template <typename Cfg>
class TPApprox {
  public:
  using real = Real<Cfg>;

  static void copyStorageToLocal(FrictionLawData<Cfg>* /*data*/,
                                 DynamicRupture::Layer& /*layerData*/,
                                 seissol::initializer::AllocationPlace /*place*/) {}

  /**
   * Use a simple copy for now, maybe use proper resampling later
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      resampleSlipRate(FrictionLawContext<Cfg, Type>& ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints<Cfg>],
                       PointValues<Cfg, Type>& resampledSlipRate) {
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      resampledSlipRate[pointIndex] = slipRateMagnitude[pointIndex];
    });
  }

  SEISSOL_HOSTDEVICE static real
      stateVariableHook(real localAccumulatedSlip, real localDc, real tpProxyExponent) {
    const real factor = (static_cast<real>(1.0) + std::fabs(localAccumulatedSlip) / localDc);
    return static_cast<real>(1.0) - std::pow(factor, -tpProxyExponent);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real strengthHook(FrictionLawContext<Cfg, Type>& /*ctx*/,
                                              real strength,
                                              real /*localSlipRate*/,
                                              real /*deltaT*/,
                                              real /*vStar*/,
                                              real /*prakashLength*/,
                                              std::uint32_t /*pointIndex*/) {
    return strength;
  }

  /**
   * See NoSpecialization::strengthHookSlope.
   */
  SEISSOL_HOSTDEVICE static real strengthHookSlope(real /*localSlipRate*/,
                                                   real /*deltaT*/,
                                                   real /*vStar*/,
                                                   real /*prakashLength*/) {
    return static_cast<real>(1.0);
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_LINEARSLIPWEAKENING_H_
