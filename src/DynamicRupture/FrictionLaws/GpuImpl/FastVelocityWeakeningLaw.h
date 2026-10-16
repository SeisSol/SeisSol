// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_

#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"

#include <cstdint>

namespace seissol::dr::friction_law::gpu {

template <typename Cfg, typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<Cfg, FastVelocityWeakeningLaw<Cfg, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;

  using RateAndStateBase<Cfg, FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data, DynamicRupture::Layer& layerData) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {
    data->srW = layerData.var<LTSRateAndStateFastVelocityWeakening::RsSrW>(
        Cfg(), seissol::initializer::AllocationPlace::Device);
  }

  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S>
  SEISSOL_DEVICE static S stateVariableAt(FrictionLawContext<Cfg>& __restrict ctx,
                                          S localSlipRate,
                                          real timeIncrement) {
    using std::expm1;
    using std::fmax;
    using std::log;
    using std::pow;
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const real localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const real localSrW = ctx.data->srW[ctx.ltsFace][ctx.pointIndex];

    const real localF0 = ctx.data->f0[ctx.ltsFace][ctx.pointIndex];
    const real localB = ctx.data->b[ctx.ltsFace][ctx.pointIndex];
    const real localMuW = ctx.data->muW[ctx.ltsFace][ctx.pointIndex];

    const S lowVelocityFriction = fmax(
        S(static_cast<real>(0)),
        S(localF0) - S(localB - localA) * log(localSlipRate / S(ctx.data->drParameters.rsSr0)));

    const S steadyStateFrictionCoefficient =
        S(localMuW) +
        (lowVelocityFriction - S(localMuW)) /
            pow(S(static_cast<real>(1.0)) + misc::power<8>(localSlipRate / S(localSrW)),
                static_cast<real>(1.0 / 8.0));

    // stated over V / (2 V_0) rather than under 2 V_0 / V: the reciprocal's derivative is
    // -2 V_0 / V^2, which leaves single precision at a slip rate where the state variable it
    // belongs to is perfectly well behaved
    const S steadyStateStateVariable =
        S(localA) *
        rs::logsinhOver(localSlipRate * S(static_cast<real>(0.5) / ctx.data->drParameters.rsSr0),
                        steadyStateFrictionCoefficient / S(localA));

    // Exact integration over the step with the slip rate held fixed. Written as psi0 + (psi_ss -
    // psi0) * (1 - exp(p)), it keeps psi_ss as its exact fixed point in any precision. The form
    // psi_ss * (1 - exp(p)) + psi0 * exp(p) rounds exp and expm1 independently, and their mismatch
    // moves the fixed point by L / (V dt) times as much, a few thousand here; in single precision
    // on a device, whose expf errs to one side, that biased the state variable of every sliding
    // point.
    const S preexp1 = -localSlipRate * S(timeIncrement / localSl0);
    const S exp1m = -expm1(preexp1);
    const S stateVarReference = S(ctx.initialVariables.stateVarReference);
    return stateVarReference + (steadyStateStateVariable - stateVarReference) * exp1m;
  }

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext<Cfg>& __restrict ctx,
                                                 real timeIncrement) {
    ctx.stateVariableBuffer =
        stateVariableAt<real>(ctx, ctx.initialVariables.localSlipRate, timeIncrement);
  }

  /// the precision the state variable of this law is stated in
  using StateScalar = real;

  /// The friction coefficient at a slip rate, with the state variable evaluated at that very slip
  /// rate. Both dependencies travel through the scalar, so a dual number comes back carrying
  /// d(mu)/dV of the composition.
  template <typename S>
  SEISSOL_DEVICE static S
      updateMuFolded(FrictionLawContext<Cfg>& __restrict ctx, S slipRate, real timeIncrement) {
    const auto stateVariable = stateVariableAt(ctx, dualCast<StateScalar>(slipRate), timeIncrement);
    const S localStateVariable = dualCast<real>(stateVariable);
    const S localA = S(ctx.data->a[ctx.ltsFace][ctx.pointIndex]);
    const S cExpLog = localStateVariable / localA;
    const S cLin = S(static_cast<real>(0.5) / ctx.data->drParameters.rsSr0);
    return localA * rs::arsinhexp(cLin * slipRate, cExpLog, rs::computeCExp(cExpLog));
  }

  /// the state variable is a closed-form function of the slip rate, so the inversion can carry
  /// it inside its own iteration instead of relaying it through a fixed point
  static constexpr bool FoldsStateVariable = true;

  struct MuDetails {
    real a{};
    real cLin{};
    real cExpLog{};
    real cExp{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext<Cfg>& __restrict ctx,
                                               real localStateVariable) {
    const real localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const real cLin = static_cast<real>(0.5) / ctx.data->drParameters.rsSr0;
    const real cExpLog = localStateVariable / localA;
    const real cExp = rs::computeCExp(cExpLog);
    return MuDetails{localA, cLin, cExpLog, cExp};
  }

  /// generic over the scalar: a dual slip rate carries the derivative out with the value
  template <typename S>
  SEISSOL_DEVICE static S updateMu(FrictionLawContext<Cfg>& /*ctx*/,
                                   S localSlipRateMagnitude,
                                   const MuDetails& details) {
    const S lx = S(details.cLin) * localSlipRateMagnitude;
    return S(details.a) * rs::arsinhexp(lx, S(details.cExpLog), S(details.cExp));
  }

  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext<Cfg>& __restrict ctx) {
    const auto localStateVariable = ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex];
    const auto toResample = ctx.stateVariableBuffer - localStateVariable;

    const auto resampledDeltaStateVar = resampleVariable(ctx, toResample);

    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] =
        localStateVariable + resampledDeltaStateVar;
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
