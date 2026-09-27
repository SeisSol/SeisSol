// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"

namespace seissol::dr::friction_law::gpu {
template <class Derived, class TPMethod>
class SlowVelocityWeakeningLaw
    : public RateAndStateBase<SlowVelocityWeakeningLaw<Derived, TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<SlowVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  // Note that we need double precision here, since single precision led to NaNs.
  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& __restrict ctx,
                                                 real timeIncrement) {
    Derived::updateStateVariable(ctx, timeIncrement);
  }

  /// the precision the state variable of this law is stated in; the relaxation rate here grows
  /// with the slip rate, so wherever the state moves the step is far above a single-precision ulp
  using StateScalar = real;

  /// The friction coefficient at a slip rate, with the state variable evaluated at that very slip
  /// rate. Both dependencies travel through the scalar, so a dual number comes back carrying
  /// d(mu)/dV of the composition.
  template <typename S>
  SEISSOL_DEVICE static S
      updateMuFolded(FrictionLawContext& __restrict ctx, S slipRate, real timeIncrement) {
    using std::log;
    const auto stateVariable =
        Derived::stateVariableAt(ctx, dualCast<StateScalar>(slipRate), timeIncrement);
    const S localStateVariable = dualCast<real>(stateVariable);
    const S localA = S(ctx.data->a[ctx.ltsFace][ctx.pointIndex]);
    const S localSl0 = S(ctx.data->sl0[ctx.ltsFace][ctx.pointIndex]);
    const S log1 = log(S(ctx.data->drParameters.rsSr0) * localStateVariable / localSl0);
    const S cExpLog = (S(ctx.data->f0[ctx.ltsFace][ctx.pointIndex]) +
                       S(ctx.data->b[ctx.ltsFace][ctx.pointIndex]) * log1) /
                      localA;
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

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext& __restrict ctx,
                                               real localStateVariable) {
    const real localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const real log1 = std::log(ctx.data->drParameters.rsSr0 * localStateVariable / localSl0);
    const real localF0 = ctx.data->f0[ctx.ltsFace][ctx.pointIndex];
    const real localB = ctx.data->b[ctx.ltsFace][ctx.pointIndex];

    const real cLin = static_cast<real>(0.5) / ctx.data->drParameters.rsSr0;
    const real cExpLog = (localF0 + localB * log1) / localA;
    const real cExp = rs::computeCExp(cExpLog);
    return MuDetails{localA, cLin, cExpLog, cExp};
  }

  /// generic over the scalar: a dual slip rate carries the derivative out with the value
  template <typename S>
  SEISSOL_DEVICE static S
      updateMu(FrictionLawContext& /*ctx*/, S localSlipRateMagnitude, const MuDetails& details) {
    const S lx = S(details.cLin) * localSlipRateMagnitude;
    return S(details.a) * rs::arsinhexp(lx, S(details.cExpLog), S(details.cExp));
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws,
   * we just copy the buffer into the member variable.
   */
  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext& __restrict ctx) {
    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] = ctx.stateVariableBuffer;
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
