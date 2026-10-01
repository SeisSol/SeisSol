// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"

namespace seissol::dr::friction_law::gpu {
template <class TPMethod>
class SevereVelocityWeakeningLaw
    : public RateAndStateBase<SevereVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<SevereVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /*
    ! friction develops as                    mu = mu_s + a V/(V+Vc) - b SV/(SV + Dc)
    ! Note the typo in eq.1 of Ampuero&Ben-Zion, 2008
    ! state variable SV develops as     dSV / dt = (V-SV) / Tc
    ! parameters: static friction mu_s, char. velocity scale Vc, charact. timescale Tc,
    ! charact. length scale Dc, direct and evolution effect coeff. a,b
    ! Notice that Dc, a and b are recycled but not equivalent to cases 3 and 4
    ! steady-state friction value is:       mu = mu_s + (a - b) V/(V+Vc)
    ! dynamic friction value (if reached) mu_d = mu_s + (a - b)
    ! Tc tunes between slip-weakening and rate-weakening behavior
    !
  */

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {}

  // Note that we need double precision here, since single precision led to NaNs.
  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S>
  SEISSOL_DEVICE static S
      stateVariableAt(FrictionLawContext& __restrict ctx, S localSlipRate, double timeIncrement) {
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];

    const S steadyStateStateVariable = localSlipRate * S(localSl0 / ctx.data->drParameters.rsSr0);

    const double preexp1 = -ctx.data->drParameters.rsSr0 * (timeIncrement / localSl0);
    const double exp1m = -std::expm1(preexp1);
    // the relaxation towards the steady state with expm1 alone, which keeps it as its exact fixed
    // point (see FastVelocityWeakeningLaw::updateStateVariable)
    const double stateVarReference = ctx.initialVariables.stateVarReference;
    return S(stateVarReference) + (steadyStateStateVariable - S(stateVarReference)) * S(exp1m);
  }

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& __restrict ctx,
                                                 double timeIncrement) {
    ctx.stateVariableBuffer = static_cast<real>(
        stateVariableAt<double>(ctx, ctx.initialVariables.localSlipRate, timeIncrement));
  }

  /*
    ! Newton-Raphson algorithm to determine the value of the slip rate.
    ! We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
    !  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i
    ).

    ! In our case we equalize the values of the traction for two equations:

    !             g =    SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))

    !             f =    (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)

    !             where mu = mu_s + a V/(V+Vc) - b SV/(SV + Vc)
  */

  /// the state variable is a closed-form function of the slip rate, so the inversion can carry
  /// it inside its own iteration instead of relaying it through a fixed point
  static constexpr bool FoldsStateVariable = true;

  /// the precision the state variable of this law is stated in. Its relaxation rate follows the
  /// reference slip rate rather than the actual one, so a step moves the state by some 2.5e-8 of
  /// itself whatever the fault is doing, which is below a single-precision ulp everywhere.
  using StateScalar = double;

  /// The friction coefficient at a slip rate, with the state variable evaluated at that very slip
  /// rate. Both dependencies travel through the scalar, so a dual number comes back carrying
  /// d(mu)/dV of the composition.
  template <typename S>
  SEISSOL_DEVICE static S
      updateMuFolded(FrictionLawContext& __restrict ctx, S slipRate, real timeIncrement) {
    const auto stateVariable =
        stateVariableAt(ctx, dualCast<StateScalar>(slipRate), static_cast<double>(timeIncrement));
    const S localStateVariable = dualCast<real>(stateVariable);
    const S localSl0 = S(ctx.data->sl0[ctx.ltsFace][ctx.pointIndex]);
    const S c = S(ctx.data->b[ctx.ltsFace][ctx.pointIndex]) * localStateVariable /
                (localStateVariable + localSl0);
    return S(ctx.data->f0[ctx.ltsFace][ctx.pointIndex]) +
           S(ctx.data->a[ctx.ltsFace][ctx.pointIndex]) * slipRate /
               (slipRate + S(ctx.data->drParameters.rsSr0)) -
           c;
  }

  struct MuDetails {
    real a{};
    real c{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext& __restrict ctx,
                                               real localStateVariable) {
    const real localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const real c = ctx.data->b[ctx.ltsFace][ctx.pointIndex] * localStateVariable /
                   (localStateVariable + localSl0);
    return MuDetails{localA, c};
  }

  /// generic over the scalar: a dual slip rate carries the derivative out with the value
  template <typename S>
  SEISSOL_DEVICE static S updateMu(FrictionLawContext& __restrict ctx,
                                   S localSlipRateMagnitude,
                                   const MuDetails& details) {
    return S(ctx.data->f0[ctx.ltsFace][ctx.pointIndex]) +
           S(details.a) * localSlipRateMagnitude /
               (localSlipRateMagnitude + S(ctx.data->drParameters.rsSr0)) -
           S(details.c);
  }

  // no resampling
  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext& __restrict ctx) {
    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] = ctx.stateVariableBuffer;
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
