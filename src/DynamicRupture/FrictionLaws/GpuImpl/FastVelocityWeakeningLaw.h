// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"

namespace seissol::dr::friction_law::gpu {

template <typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {}

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {
    data->srW = layerData.var<LTSRateAndStateFastVelocityWeakening::RsSrW>(
        seissol::initializer::AllocationPlace::Device);
  }

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& ctx, real timeIncrement) {
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const real localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const real localSrW = ctx.data->srW[ctx.ltsFace][ctx.pointIndex];
    const real localSlipRate = ctx.initialVariables.localSlipRate;

    const auto localF0 = ctx.data->f0[ctx.ltsFace][ctx.pointIndex];
    const auto localB = ctx.data->b[ctx.ltsFace][ctx.pointIndex];
    const auto localMuW = ctx.data->muW[ctx.ltsFace][ctx.pointIndex];

    const real lowVelocityFriction =
        localF0 - (localB - localA) * std::log(localSlipRate / ctx.data->drParameters.rsSr0);

    const real steadyStateFrictionCoefficient =
        localMuW + (lowVelocityFriction - localMuW) /
                       std::pow(1.0 + std::pow(localSlipRate / localSrW, 8), 1.0 / 8.0);

    const real steadyStateStateVariable =
        localA * rs::logsinh(ctx.data->drParameters.rsSr0 / localSlipRate * 2,
                             steadyStateFrictionCoefficient / localA);

    const double preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const double exp1v = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    const real localStateVariable =
        steadyStateStateVariable * exp1m + exp1v * ctx.initialVariables.stateVarReference;

    ctx.stateVariableBuffer = localStateVariable;
  }

  struct MuDetails {
    real a{};
    real cLin{};
    real cExpLog{};
    real cExp{};
    real acLin{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext& ctx, double localStateVariable) {
    const real localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const real cLin = 0.5 / ctx.data->drParameters.rsSr0;
    const real cExpLog = localStateVariable / localA;
    const real cExp = rs::computeCExp(cExpLog);
    const real acLin = localA * cLin;
    return MuDetails{localA, cLin, cExpLog, cExp, acLin};
  }

  SEISSOL_DEVICE static real
      updateMu(FrictionLawContext& ctx, real localSlipRateMagnitude, const MuDetails& details) {
    const real lx = details.cLin * localSlipRateMagnitude;
    return details.a * rs::arsinhexp(lx, details.cExpLog, details.cExp);
  }

  SEISSOL_DEVICE static real updateMuDerivative(FrictionLawContext& ctx,
                                                real localSlipRateMagnitude,
                                                const MuDetails& details) {
    const real lx = details.cLin * localSlipRateMagnitude;
    return details.acLin * rs::arsinhexpDerivative(lx, details.cExpLog, details.cExp);
  }

  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext& ctx) {
    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    static_assert(Dim0 == misc::NumPaddedPointsSingleSim);
    static_assert(Dim0 >= Dim1);

    const auto localStateVariable = ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex];
    ctx.sharedMemory[ctx.pointIndex] = ctx.stateVariableBuffer - localStateVariable;
    deviceBarrier(ctx);

    const auto simPointIndex = ctx.pointIndex / multisim::NumSimulations;
    const auto simId = ctx.pointIndex % multisim::NumSimulations;

    real resampledDeltaStateVar{0.0};
    for (size_t i{0}; i < Dim1; ++i) {
      if constexpr (multisim::MultisimEnabled) {
        resampledDeltaStateVar += ctx.args->resampleMatrix[simPointIndex * Dim1 + i] *
                                  ctx.sharedMemory[i * multisim::NumSimulations + simId];
      } else {
        resampledDeltaStateVar +=
            ctx.args->resampleMatrix[simPointIndex + i * Dim0] * ctx.sharedMemory[i];
      }
    }

    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] =
        localStateVariable + resampledDeltaStateVar;
  }

  SEISSOL_DEVICE static void executeIfNotConverged() {}
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
