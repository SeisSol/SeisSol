// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::friction_law::gpu {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <typename Derived>
class LinearSlipWeakeningBase : public BaseFrictionSolver<LinearSlipWeakeningBase<Derived>> {
  public:
  explicit LinearSlipWeakeningBase(const FrictionLawParameters& drParameters)
      : BaseFrictionSolver<LinearSlipWeakeningBase<Derived>>(drParameters) {};

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {
    Derived::copySpecificStorageDataToLocal(data, layerData);
  }

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& __restrict ctx,
                                                   uint32_t timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    Derived::calcStrengthHook(ctx, timeIndex);
    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    calcSlipRateAndTraction(ctx, timeIndex);
    Derived::calcStateVariableHook(ctx, timeIndex);
    frictionFunctionHook(ctx);
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  SEISSOL_DEVICE static void calcSlipRateAndTraction(FrictionLawContext& __restrict ctx,
                                                     uint32_t timeIndex) {
    const auto deltaT{ctx.args->deltaT[timeIndex]};

    auto& faultStresses = ctx.faultStresses;
    auto& tractionResults = ctx.tractionResults;
    auto& strength = ctx.strengthBuffer;

    // calculate absolute value of stress in Y and Z direction
    const real totalStress1 =
        ctx.data->initialStressInFaultCS[ctx.ltsFace][3][ctx.pointIndex] + faultStresses.traction1;
    const real totalStress2 =
        ctx.data->initialStressInFaultCS[ctx.ltsFace][5][ctx.pointIndex] + faultStresses.traction2;
    const real absoluteShearStress = misc::magnitude(totalStress1, totalStress2);

    const auto [eta, invEta] = common::projectEta(ctx.data->impAndEta[ctx.ltsFace],
                                                  ctx.data->impedanceMatrices[ctx.ltsFace],
                                                  totalStress1,
                                                  totalStress2,
                                                  absoluteShearStress);

    // the direction along which the slip rate is decomposed further down; scaled such that
    // dividing by `divisor` yields the unit slip direction
    real dirStress1 = totalStress1;
    real dirStress2 = totalStress2;
    real etaEff = eta;
    real slipRateMagnitude{};

    if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
      // see the CPU implementation for the derivation: sweeping n -> V -> n twice resolves both
      // the non-collinearity of slip and trial traction and the normal/shear coupling
      const real invAbsolute = (absoluteShearStress > 0)
                                   ? static_cast<real>(1.0) / absoluteShearStress
                                   : static_cast<real>(0.0);
      real n1 = totalStress1 * invAbsolute;
      real n2 = totalStress2 * invAbsolute;
      real projectedStress = absoluteShearStress;
      real localEta = eta;
      real localEtaNormal = common::projectEtaNormal(ctx.data->impAndEta[ctx.ltsFace],
                                                     ctx.data->impedanceMatrices[ctx.ltsFace],
                                                     totalStress1,
                                                     totalStress2,
                                                     absoluteShearStress);

      constexpr std::uint32_t DirectionSweeps = 2;
      for (std::uint32_t sweep = 0; sweep < DirectionSweeps; ++sweep) {
        etaEff = localEta + ctx.strengthSlopeBuffer * localEtaNormal;
        etaEff = (etaEff > 0) ? etaEff : localEta;
        slipRateMagnitude = std::max(static_cast<real>(0.0), (projectedStress - strength) / etaEff);

        if (sweep + 1 == DirectionSweeps) {
          break;
        }

        const real localStrength = projectedStress - slipRateMagnitude * localEta;
        const auto [d1, d2] = common::updateSlipDirection(ctx.data->impAndEta[ctx.ltsFace],
                                                          ctx.data->impedanceMatrices[ctx.ltsFace],
                                                          localStrength,
                                                          slipRateMagnitude,
                                                          totalStress1,
                                                          totalStress2,
                                                          absoluteShearStress);
        n1 = d1;
        n2 = d2;
        projectedStress = n1 * totalStress1 + n2 * totalStress2;
        const auto [e, unusedInv] = common::projectEta(ctx.data->impAndEta[ctx.ltsFace],
                                                       ctx.data->impedanceMatrices[ctx.ltsFace],
                                                       n1,
                                                       n2,
                                                       static_cast<real>(1.0));
        localEta = e;
        localEtaNormal = common::projectEtaNormal(ctx.data->impAndEta[ctx.ltsFace],
                                                  ctx.data->impedanceMatrices[ctx.ltsFace],
                                                  n1,
                                                  n2,
                                                  static_cast<real>(1.0));
      }

      dirStress1 = n1 * projectedStress;
      dirStress2 = n2 * projectedStress;
    } else {
      slipRateMagnitude =
          std::max(static_cast<real>(0.0), (absoluteShearStress - strength) * invEta);
    }

    // calculate slip rates
    ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] = slipRateMagnitude;
    const auto divisor = strength + etaEff * slipRateMagnitude;
    ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] = slipRateMagnitude * dirStress1 / divisor;
    ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] = slipRateMagnitude * dirStress2 / divisor;

    const auto [tU1, tU2] = common::matmulEta(ctx.data->impAndEta[ctx.ltsFace],
                                              ctx.data->impedanceMatrices[ctx.ltsFace],
                                              ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex],
                                              ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex]);

    const auto tUN = common::matmulEtaNormal(ctx.data->impAndEta[ctx.ltsFace],
                                             ctx.data->impedanceMatrices[ctx.ltsFace],
                                             ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex],
                                             ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex]);

    // calculate traction
    // the normal stress written here is the *dynamic* normal traction, using the slip rate just
    // solved for
    tractionResults.normalStress = faultStresses.normalStress - tUN;
    tractionResults.traction1 = faultStresses.traction1 - tU1;
    tractionResults.traction2 = faultStresses.traction2 - tU2;
    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] = tractionResults.traction1;
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] = tractionResults.traction2;
    // update directional slip
    ctx.data->slip1[ctx.ltsFace][ctx.pointIndex] +=
        ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] * deltaT;
    ctx.data->slip2[ctx.ltsFace][ctx.pointIndex] +=
        ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] * deltaT;
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  SEISSOL_DEVICE static void frictionFunctionHook(FrictionLawContext& __restrict ctx) {
    auto& stateVariable = ctx.stateVariableBuffer;
    ctx.data->mu[ctx.ltsFace][ctx.pointIndex] =
        ctx.data->muS[ctx.ltsFace][ctx.pointIndex] -
        (ctx.data->muS[ctx.ltsFace][ctx.pointIndex] - ctx.data->muD[ctx.ltsFace][ctx.pointIndex]) *
            stateVariable;
    // instantaneous healing
    if ((ctx.data->peakSlipRate[ctx.ltsFace][ctx.pointIndex] >
         ctx.data->drParameters.healingThreshold) &&
        (ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] <
         ctx.data->drParameters.healingThreshold)) {
      ctx.data->mu[ctx.ltsFace][ctx.pointIndex] = ctx.data->muS[ctx.ltsFace][ctx.pointIndex];
      stateVariable = 0.0;
    }
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& __restrict ctx,
                                                     real time) {
    if (ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] &&
        std::fabs(ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex]) >=
            ctx.data->dC[ctx.ltsFace][ctx.pointIndex]) {
      ctx.data->dynStressTime[ctx.ltsFace][ctx.pointIndex] = time;
      ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] = false;
    }
  }

  SEISSOL_DEVICE static void preHook(FrictionLawContext& __restrict ctx) {}
  SEISSOL_DEVICE static void postHook(FrictionLawContext& __restrict ctx) {}

  protected:
  static constexpr real U0 = 10e-14;
};

template <class SpecializationT>
class LinearSlipWeakeningLaw
    : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>> {
  public:
  explicit LinearSlipWeakeningLaw(const FrictionLawParameters& drParameters)
      : LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>>(drParameters),
        specialization_(drParameters) {};

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {
    data->dC =
        layerData.var<LTSLinearSlipWeakening::DC>(seissol::initializer::AllocationPlace::Device);
    data->muS =
        layerData.var<LTSLinearSlipWeakening::MuS>(seissol::initializer::AllocationPlace::Device);
    data->muD =
        layerData.var<LTSLinearSlipWeakening::MuD>(seissol::initializer::AllocationPlace::Device);
    data->cohesion = layerData.var<LTSLinearSlipWeakening::Cohesion>(
        seissol::initializer::AllocationPlace::Device);
    data->forcedRuptureTime = layerData.var<LTSLinearSlipWeakening::ForcedRuptureTime>(
        seissol::initializer::AllocationPlace::Device);
    SpecializationT::copyStorageToLocal(data, layerData);
  }

  SEISSOL_DEVICE static void calcStrengthHook(FrictionLawContext& __restrict ctx,
                                              uint32_t timeIndex) {

    const auto deltaT{ctx.args->deltaT[timeIndex]};

    const auto vStar{ctx.data->drParameters.vStar};
    const auto prakashLength{ctx.data->drParameters.prakashLength};

    auto& strength = ctx.strengthBuffer;

    // The anisotropic normal/shear coupling is deliberately not applied here: the strength is
    // affine in the normal stress, so it is handled exactly through the divisor in
    // calcSlipRateAndTraction, using the slope filled in below.
    const real totalNormalStress =
        ctx.data->initialStressInFaultCS[ctx.ltsFace][0][ctx.pointIndex] +
        ctx.faultStresses.normalStress + ctx.data->initialPressure[ctx.ltsFace][ctx.pointIndex] +
        ctx.faultStresses.fluidPressure;
    strength = -ctx.data->cohesion[ctx.ltsFace][ctx.pointIndex] -
               ctx.data->mu[ctx.ltsFace][ctx.pointIndex] *
                   std::min(totalNormalStress, static_cast<real>(0.0));

    strength =
        SpecializationT::strengthHook(ctx,
                                      strength,
                                      ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex],
                                      deltaT,
                                      vStar,
                                      prakashLength);

    if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
      // d(strength) / d(-sigma_eff). Zero while the normal stress is clamped; the clamp is
      // evaluated at the uncorrected normal stress, which is second order in the coupling.
      ctx.strengthSlopeBuffer = (totalNormalStress < 0 ? ctx.data->mu[ctx.ltsFace][ctx.pointIndex]
                                                       : static_cast<real>(0.0)) *
                                SpecializationT::strengthHookSlope(
                                    ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex],
                                    deltaT,
                                    vStar,
                                    prakashLength);
    }
  }

  SEISSOL_DEVICE static void calcStateVariableHook(FrictionLawContext& __restrict ctx,
                                                   uint32_t timeIndex) {
    const auto t0{ctx.data->drParameters.t0[0]};
    const auto tpProxyExponent{ctx.data->drParameters.tpProxyExponent};

    real tn = ctx.args->fullUpdateTime;
    for (uint32_t i = 0; i <= timeIndex; ++i) {
      tn += ctx.args->deltaT[i];
    }

    const real resampledSlipRate =
        SpecializationT::resampleSlipRate(ctx, ctx.data->slipRateMagnitude[ctx.ltsFace]);

    // integrate slip rate to get slip = state variable
    ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex] +=
        resampledSlipRate * ctx.args->deltaT[timeIndex];

    // Actually slip is already the stateVariable for this FL, but to simplify the next
    // equations we divide it here by the critical distance.
    const real localStateVariable = SpecializationT::stateVariableHook(
        ctx,
        ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex],
        ctx.data->dC[ctx.ltsFace][ctx.pointIndex],
        tpProxyExponent);

    real f2 = 0.0;
    if (t0 == 0) {
      f2 = static_cast<real>(tn >= ctx.data->forcedRuptureTime[ctx.ltsFace][ctx.pointIndex]);
    } else {
      f2 = misc::clamp((tn - ctx.data->forcedRuptureTime[ctx.ltsFace][ctx.pointIndex]) / t0,
                       static_cast<real>(0.0),
                       static_cast<real>(1.0));
    }
    ctx.stateVariableBuffer = std::max(localStateVariable, f2);
  }

  protected:
  SpecializationT specialization_;
};

class NoSpecialization {
  public:
  explicit NoSpecialization(const FrictionLawParameters& parameters) {};

  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {}

  SEISSOL_DEVICE static real
      resampleSlipRate(FrictionLawContext& __restrict ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
    return resampleVariable(ctx, slipRateMagnitude[ctx.pointIndex]);
  };

  SEISSOL_DEVICE static real stateVariableHook(FrictionLawContext& /*ctx*/,
                                               real localAccumulatedSlip,
                                               real localDc,
                                               real /*tpProxyExponent*/) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  SEISSOL_DEVICE static real strengthHook(FrictionLawContext& /*ctx*/,
                                          real strength,
                                          real /*localSlipRate*/,
                                          real /*deltaT*/,
                                          real /*vStar*/,
                                          real /*prakashLength*/) {
    return strength;
  };

  /**
   * d(strengthHook output) / d(its faultStrength argument). Only needed for the anisotropic
   * normal/shear coupling. MUST be free of side effects and evaluated with the same arguments as
   * the corresponding strengthHook call.
   */
  SEISSOL_DEVICE static real strengthHookSlope(real /*localSlipRate*/,
                                               real /*deltaT*/,
                                               real /*vStar*/,
                                               real /*prakashLength*/) {
    return static_cast<real>(1.0);
  };
};

class BiMaterialFault {
  public:
  explicit BiMaterialFault(const FrictionLawParameters& parameters) {};

  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {
    data->regularizedStrength =
        layerData.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(
            seissol::initializer::AllocationPlace::Device);
  }

  SEISSOL_DEVICE static real
      resampleSlipRate(FrictionLawContext& __restrict ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
    return slipRateMagnitude[ctx.pointIndex];
  };

  SEISSOL_DEVICE static real stateVariableHook(FrictionLawContext& /*ctx*/,
                                               real localAccumulatedSlip,
                                               real localDc,
                                               real /*tpProxyExponent*/) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  SEISSOL_DEVICE static real strengthHook(FrictionLawContext& __restrict ctx,
                                          real faultStrength,
                                          real localSlipRate,
                                          real deltaT,
                                          real vStar,
                                          real prakashLength) {
    const auto expval =
        -(std::max(static_cast<real>(0.0), localSlipRate) + vStar) * deltaT / prakashLength;
    const real expterm = std::exp(expval);
    const real exp1mterm = -std::expm1(expval);

    const real newStrength = ctx.data->regularizedStrength[ctx.ltsFace][ctx.pointIndex] * expterm +
                             faultStrength * exp1mterm;

    ctx.data->regularizedStrength[ctx.ltsFace][ctx.pointIndex] = newStrength;
    return newStrength;
  };

  /**
   * See NoSpecialization::strengthHookSlope. The Prakash-Clifton regularisation low-passes the
   * strength, so only the fraction exp1mterm of a change in faultStrength arrives instantaneously;
   * regularizedStrength carries the previous step and drops out of the derivative.
   */
  SEISSOL_DEVICE static real
      strengthHookSlope(real localSlipRate, real deltaT, real vStar, real prakashLength) {
    const auto expval =
        -(std::max(static_cast<real>(0.0), localSlipRate) + vStar) * deltaT / prakashLength;
    return -std::expm1(expval);
  };
};

class TPApprox {
  public:
  explicit TPApprox(const FrictionLawParameters& parameters) {};

  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {}

  SEISSOL_DEVICE static real
      resampleSlipRate(FrictionLawContext& __restrict ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
    return slipRateMagnitude[ctx.pointIndex];
  };

  SEISSOL_DEVICE static real stateVariableHook(FrictionLawContext& /*ctx*/,
                                               real localAccumulatedSlip,
                                               real localDc,
                                               real tpProxyExponent) {
    const real factor = (static_cast<real>(1.0) + std::fabs(localAccumulatedSlip) / localDc);
    return static_cast<real>(1.0) - std::pow(factor, -tpProxyExponent);
  };

  SEISSOL_DEVICE static real strengthHook(FrictionLawContext& /*ctx*/,
                                          real strength,
                                          real /*localSlipRate*/,
                                          real /*deltaT*/,
                                          real /*vStar*/,
                                          real /*prakashLength*/) {
    return strength;
  };

  /**
   * d(strengthHook output) / d(its faultStrength argument). Only needed for the anisotropic
   * normal/shear coupling. MUST be free of side effects and evaluated with the same arguments as
   * the corresponding strengthHook call.
   */
  SEISSOL_DEVICE static real strengthHookSlope(real /*localSlipRate*/,
                                               real /*deltaT*/,
                                               real /*vStar*/,
                                               real /*prakashLength*/) {
    return static_cast<real>(1.0);
  };
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_
