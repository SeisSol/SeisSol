// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SEVEREVELOCITYWEAKENINGLAW_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Physics/Dual.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/RateAndState.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cmath>
#include <cstdint>

namespace seissol::dr::friction_law {

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
template <typename Cfg, typename TPMethod>
class SevereVelocityWeakeningLaw
    : public RateAndState<Cfg, SevereVelocityWeakeningLaw<Cfg, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;

  static void copyLawStorageToLocal(FrictionLawData<Cfg>* /*data*/,
                                    DynamicRupture::Layer& /*layerData*/,
                                    seissol::initializer::AllocationPlace /*place*/) {}

  /// the precision the state variable of this law is stated in. Unlike the other laws the
  /// relaxation rate here is -expm1(-V0 dt / L): it follows the reference slip rate rather than
  /// the actual one, so a step moves the state by some 2.5e-8 of itself whatever the fault is
  /// doing, which is below a single-precision ulp at every slip rate.
  using StateScalar = double;

  /// the state variable is a closed-form function of the slip rate, so the inversion can carry it
  /// inside its own iteration instead of relaying it through a fixed point
  static constexpr bool FoldsStateVariable = true;

  // Note that we need double precision here, since single precision led to NaNs.
  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S stateVariableAt(const FrictionLawContext<Cfg, Type>& ctx,
                                              std::uint32_t pointIndex,
                                              double stateVarReference,
                                              double timeIncrement,
                                              S localSlipRate) {
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][pointIndex];

    const S steadyStateStateVariable = localSlipRate * S(localSl0 / ctx.data->drParameters.rsSr0);

    const double preexp1 = -ctx.data->drParameters.rsSr0 * (timeIncrement / localSl0);
    const double exp1m = -std::expm1(preexp1);
    // the relaxation towards the steady state with expm1 alone, which keeps it as its exact fixed
    // point (see FastVelocityWeakeningLaw::stateVariableAt)
    return S(stateVarReference) + (steadyStateStateVariable - S(stateVarReference)) * S(exp1m);
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

  /// The friction coefficient at a slip rate, with the state variable evaluated at that very slip
  /// rate. Both dependencies travel through the scalar, so a dual number comes back carrying
  /// d(mu)/dV of the composition.
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S updateMuFolded(const FrictionLawContext<Cfg, Type>& ctx,
                                             std::uint32_t pointIndex,
                                             S slipRate,
                                             real stateVarReference,
                                             real timeIncrement) {
    const auto ltsFace = ctx.ltsFace;
    const auto stateVariable = stateVariableAt(ctx,
                                               pointIndex,
                                               static_cast<double>(stateVarReference),
                                               static_cast<double>(timeIncrement),
                                               dualCast<StateScalar>(slipRate));

    const S localStateVariable = dualCast<real>(stateVariable);
    const S localSl0 = S(ctx.data->sl0[ltsFace][pointIndex]);
    const S c =
        S(ctx.data->b[ltsFace][pointIndex]) * localStateVariable / (localStateVariable + localSl0);
    return S(ctx.data->f0[ltsFace][pointIndex]) +
           S(ctx.data->a[ltsFace][pointIndex]) * slipRate /
               (slipRate + S(ctx.data->drParameters.rsSr0)) -
           c;
  }

  template <common::RangeType Type>
  struct MuDetails {
    PointValues<Cfg, Type> a{};
    PointValues<Cfg, Type> c{};
    PointValues<Cfg, Type> f0{};
  };

  /// the coefficients of mu precomputed from the state variable (ctx.stateVariableBuffer) of
  /// the points of the range
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static MuDetails<Type> getMuDetails(const FrictionLawContext<Cfg, Type>& ctx) {
    MuDetails<Type> details{};
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      const real localStateVariable = ctx.stateVariableBuffer[pointIndex];
      const auto ltsFace = ctx.ltsFace;
      const real localA = ctx.data->a[ltsFace][pointIndex];
      const real localSl0 = ctx.data->sl0[ltsFace][pointIndex];
      const real c =
          ctx.data->b[ltsFace][pointIndex] * localStateVariable / (localStateVariable + localSl0);

      details.a[pointIndex] = localA;
      details.c[pointIndex] = c;
      details.f0[pointIndex] = ctx.data->f0[ltsFace][pointIndex];
    });
    return details;
  }

  /// generic over the scalar: a dual slip rate carries the derivative out with the value
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S updateMu(const FrictionLawContext<Cfg, Type>& ctx,
                                       std::uint32_t pointIndex,
                                       S localSlipRateMagnitude,
                                       const MuDetails<Type>& details) {
    return S(details.f0[pointIndex]) +
           S(details.a[pointIndex]) * localSlipRateMagnitude /
               (localSlipRateMagnitude + S(ctx.data->drParameters.rsSr0)) -
           S(details.c[pointIndex]);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws, we just copy the buffer into the
   * member variable.
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void resampleStateVar(FrictionLawContext<Cfg, Type>& ctx) {
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      ctx.data->stateVariable[ctx.ltsFace][pointIndex] = ctx.stateVariableBuffer[pointIndex];
    });
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SEVEREVELOCITYWEAKENINGLAW_H_
