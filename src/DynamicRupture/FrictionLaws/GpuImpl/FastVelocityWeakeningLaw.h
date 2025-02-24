// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

namespace seissol::dr::friction_law::gpu {

template <typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {}

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {
    using SelfInitializerType = seissol::initializer::LTSRateAndStateFastVelocityWeakening;
    const auto* concreteLts = dynamic_cast<const SelfInitializerType*>(dynRup);
    data->srW = layerData.var(concreteLts->rsSrW, seissol::initializer::AllocationPlace::Device);
  }

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& ctx, double timeIncrement) {
    const double muW{ctx.data->drParameters.muW};

    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const double localSrW = ctx.data->srW[ctx.ltsFace][ctx.pointIndex];
    const double localSlipRate = ctx.initialVariables.localSlipRate;

    const double lowVelocityFriction =
        ctx.data->drParameters.rsF0 - (ctx.data->drParameters.rsB - localA) *
                                          std::log(localSlipRate / ctx.data->drParameters.rsSr0);

    const double steadyStateFrictionCoefficient =
        muW + (lowVelocityFriction - muW) /
                  std::pow(1.0 + std::pow(localSlipRate / localSrW, 8), 1.0 / 8.0);

    const double steadyStateStateVariable =
        localA * std::log(ctx.data->drParameters.rsSr0 / localSlipRate * 2 *
                          std::sinh(steadyStateFrictionCoefficient / localA));

    const double preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const double exp1 = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    const double localStateVariable =
        steadyStateStateVariable * exp1m + exp1 * ctx.initialVariables.stateVarReference;

    ctx.stateVariableBuffer = localStateVariable;
  }

  struct MuDetails {
    double a{};
    double c{};
    double ac{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext& ctx, double localStateVariable) {
    const double localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const double c = 0.5 / ctx.data->drParameters.rsSr0 * std::exp(localStateVariable / localA);
    return MuDetails{localA, c, localA * c};
  }

  SEISSOL_DEVICE static double
      updateMu(FrictionLawContext& ctx, double localSlipRateMagnitude, MuDetails& details) {
    const double x = details.c * localSlipRateMagnitude;
    return details.a * std::asinh(x);
  }

  SEISSOL_DEVICE static double updateMuDerivative(FrictionLawContext& ctx,
                                                  double localSlipRateMagnitude,
                                                  MuDetails& details) {
    const double x = details.c * localSlipRateMagnitude;
    return details.ac / std::sqrt(std::pow(x, 2) + 1.0);
  }

  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext& ctx) {
    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    static_assert(Dim0 == misc::NumPaddedPoints);
    static_assert(Dim0 >= Dim1);

    const auto localStateVariable = ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex];
    ctx.sharedMemory[ctx.pointIndex] = ctx.stateVariableBuffer - localStateVariable;
    deviceBarrier(ctx);

    real resampledDeltaStateVar{0.0};
    for (size_t i{0}; i < Dim1; ++i) {
      resampledDeltaStateVar += ctx.resampleMatrix[ctx.pointIndex + i * Dim0] * ctx.sharedMemory[i];
    }

    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] =
        localStateVariable + resampledDeltaStateVar;
  }

  SEISSOL_DEVICE static void executeIfNotConverged() {}
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
