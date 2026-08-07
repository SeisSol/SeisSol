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
      const auto absoluteShearStress = ctx.initialVariables.absoluteShearTraction;
      const auto localSlipRateMagnitude = ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];

      real slipRateTest{0};
      real exportMu{0};

      const bool hasConvergedLocal =
          RateAndStateBase::invertSlipRateIterative(ctx,
                                                    slipRateTest,
                                                    localStateVariable,
                                                    normalStress,
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
    // `divisor` yields the unit slip direction. For isotropy slipDirection is the normalised trial
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

  SEISSOL_DEVICE static bool invertSlipRateIterative(FrictionLawContext& __restrict ctx,
                                                     real& slipRateTest,
                                                     real localStateVariable,
                                                     real normalStress,
                                                     real absoluteShearStress,
                                                     real slipRateMagnitude,
                                                     real invEtaS,
                                                     real& exportMu) {

    // Note that we need double precision here, since single precision led to NaNs.
    real muF{0.0};
    real dMuF{0.0};
    real g{0.0};
    real dG{0.0};
    slipRateTest = slipRateMagnitude;

    const auto details = Derived::getMuDetails(ctx, localStateVariable);

    for (uint32_t i = 0; i < ctx.data->drParameters.rsMaxNumberSlipRateUpdates; i++) {
      muF = Derived::updateMu(ctx, slipRateTest, details);

      g = -invEtaS * (std::abs(normalStress) * muF - absoluteShearStress) - slipRateTest;

      const bool converged = std::abs(g) < ctx.data->drParameters.rsSlipRateTolerance;

      if (converged) {
        // we've reached the fixed point
        // NOTE: in doubt, a fixed-point mu can be recovered from slipRateTest at this point.
        // just invert -invEtaS * (std::abs(normalStress) * muF - absoluteShearStress) ==
        // slipRateTest for muF in that case.
        exportMu = muF;
        return true;
      }

      dMuF = Derived::updateMuDerivative(ctx, slipRateTest, details);
      dG = -invEtaS * (std::abs(normalStress) * dMuF) - static_cast<real>(1.0);
      slipRateTest = std::max(friction_law::rs::almostZero(), slipRateTest - (g / dG));
    }
    return false;
  }

  /**
   * Effective normal stress, including the anisotropic normal/shear coupling.
   *
   * With the 3x3 impedance eta, shear slip changes the fault-normal traction:
   *   sigma(V) = sigma_stick - V * etaNormal,   etaNormal = (eta * n)_n
   * For every isotropic material etaNormal is zero and this reduces to the previous formula.
   *
   * The slip rate is taken from the previous outer fixed-point iteration (or, on entry, from the
   * previous time step), so the Newton solver itself still sees a frozen sigma.
   */
  SEISSOL_DEVICE static void updateNormalStress(FrictionLawContext& __restrict ctx) {
    ctx.initialVariables.normalStress =
        std::min(static_cast<real>(0.0),
                 ctx.faultStresses.normalStress +
                     ctx.data->initialStressInFaultCS[ctx.ltsFace][0][ctx.pointIndex] +
                     ctx.faultStresses.fluidPressure +
                     ctx.data->initialPressure[ctx.ltsFace][ctx.pointIndex] -
                     TPMethod::getFluidPressure(ctx) -
                     ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] *
                         ctx.initialVariables.etaNormal);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_RATEANDSTATE_H_
