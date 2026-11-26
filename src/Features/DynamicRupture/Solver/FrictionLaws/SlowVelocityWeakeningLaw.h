// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SLOWVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SLOWVELOCITYWEAKENINGLAW_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Physics/Dual.h"
#include "Features/DynamicRupture/Physics/RateAndStateCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/RateAndState.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cmath>
#include <cstdint>

namespace seissol::dr::friction_law {

/**
 * The rate-and-state laws with the logarithmic friction coefficient of the slow velocity weakening,
 * once for the host and the device; `Derived` (AgingLaw, SlipLaw) supplies the evolution of the
 * state variable (stateVariableAt).
 */
template <typename Cfg, typename Derived, typename TPMethod>
class SlowVelocityWeakeningLaw
    : public RateAndState<Cfg, SlowVelocityWeakeningLaw<Cfg, Derived, TPMethod>, TPMethod> {
  private:
  // only the CRTP subclass `Derived` constructs this base
  SlowVelocityWeakeningLaw() = default;
  friend Derived;

  public:
  using real = Real<Cfg>;

  static void copyLawStorageToLocal(FrictionLawData<Cfg>* /*data*/,
                                    DynamicRupture::Layer& /*layerData*/,
                                    seissol::initializer::AllocationPlace /*place*/) {}

  /// the precision the state variable of this law is stated in. The relaxation rate here is
  /// -expm1(-V dt / L) and grows with the slip rate, so wherever the state actually moves the step
  /// is far above a single-precision ulp: some 8e-3 of the state at a metre per second.
  using StateScalar = real;

  /// the state variable is a closed-form function of the slip rate, so the inversion can carry it
  /// inside its own iteration instead of relaying it through a fixed point
  static constexpr bool FoldsStateVariable = true;

  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S stateVariableAt(const FrictionLawContext<Cfg, Type>& ctx,
                                              std::uint32_t pointIndex,
                                              real stateVarReference,
                                              real timeIncrement,
                                              S localSlipRate) {
    return Derived::stateVariableAt(
        ctx, pointIndex, stateVarReference, timeIncrement, localSlipRate);
  }

  /// The friction coefficient at a slip rate, with the state variable evaluated at that very slip
  /// rate. Both dependencies travel through the scalar, so a dual number comes back carrying
  /// d(mu)/dV of the composition.
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S updateMuFolded(const FrictionLawContext<Cfg, Type>& ctx,
                                             std::uint32_t pointIndex,
                                             S slipRate,
                                             real stateVarReference,
                                             real timeIncrement) {
    using std::log;
    const auto ltsFace = ctx.ltsFace;
    const auto stateVariable = Derived::stateVariableAt(
        ctx, pointIndex, stateVarReference, timeIncrement, dualCast<StateScalar>(slipRate));

    const S localStateVariable = dualCast<real>(stateVariable);
    const S localA = S(ctx.data->a[ltsFace][pointIndex]);
    const S localSl0 = S(ctx.data->sl0[ltsFace][pointIndex]);
    const S log1 = log(S(ctx.data->drParameters.rsSr0) * localStateVariable / localSl0);
    const S cExpLog =
        (S(ctx.data->f0[ltsFace][pointIndex]) + S(ctx.data->b[ltsFace][pointIndex]) * log1) /
        localA;
    const S cLin = S(static_cast<real>(0.5) / ctx.data->drParameters.rsSr0);
    return localA * rs::arsinhexp(cLin * slipRate, cExpLog, rs::computeCExp(cExpLog));
  }

  template <common::RangeType Type>
  struct MuDetails {
    PointValues<Cfg, Type> a{};
    PointValues<Cfg, Type> cLin{};
    PointValues<Cfg, Type> cExpLog{};
    PointValues<Cfg, Type> cExp{};
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
      const real log1 = std::log(ctx.data->drParameters.rsSr0 * localStateVariable / localSl0);
      const real localF0 = ctx.data->f0[ltsFace][pointIndex];
      const real localB = ctx.data->b[ltsFace][pointIndex];

      const real cLin = static_cast<real>(0.5) / ctx.data->drParameters.rsSr0;
      const real cExpLog = (localF0 + localB * log1) / localA;
      const real cExp = rs::computeCExp(cExpLog);

      details.a[pointIndex] = localA;
      details.cLin[pointIndex] = cLin;
      details.cExpLog[pointIndex] = cExpLog;
      details.cExp[pointIndex] = cExp;
    });
    return details;
  }

  /**
   * Computes the friction coefficient from the state variable and slip rate
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Psi
   * / L)}{a} \right)\right).\f]
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Psi \f$
   * @return \f$ \mu \f$
   */
  /// generic over the scalar: a dual slip rate carries the derivative out with the value
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S updateMu(const FrictionLawContext<Cfg, Type>& /*ctx*/,
                                       std::uint32_t pointIndex,
                                       S localSlipRateMagnitude,
                                       const MuDetails<Type>& details) {
    const S lx = S(details.cLin[pointIndex]) * localSlipRateMagnitude;
    return S(details.a[pointIndex]) *
           rs::arsinhexp(lx, S(details.cExpLog[pointIndex]), S(details.cExp[pointIndex]));
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

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SLOWVELOCITYWEAKENINGLAW_H_
