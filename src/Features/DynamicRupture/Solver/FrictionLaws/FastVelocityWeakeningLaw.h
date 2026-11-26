// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FASTVELOCITYWEAKENINGLAW_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Physics/Dual.h"
#include "Features/DynamicRupture/Physics/RateAndStateCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/RateAndState.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cassert>
#include <cmath>
#include <cstdint>

namespace seissol::dr::friction_law {

template <typename Cfg, typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndState<Cfg, FastVelocityWeakeningLaw<Cfg, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;

  /**
   * Copies all parameters from the DynamicRupture LTS to the local attributes
   */
  static void copyLawStorageToLocal(FrictionLawData<Cfg>* data,
                                    DynamicRupture::Layer& layerData,
                                    seissol::initializer::AllocationPlace place) {
    data->srW = layerData.var<LTSRateAndStateFastVelocityWeakening::RsSrW>(Cfg(), place);
  }

  /// the precision the state variable of this law is stated in
  using StateScalar = real;

  /// the state variable is a closed-form function of the slip rate, so the inversion can carry it
  /// inside its own iteration instead of relaying it through a fixed point
  static constexpr bool FoldsStateVariable = true;

  /**
   * Integrates the state variable ODE in time
   * \f[\frac{\partial \Psi}{\partial t} = - \frac{V}{L}\left(\Psi - \Psi_{ss}(V) \right)\f]
   * with steady state variable \f$\Psi_{ss}\f$.
   * Assume \f$V\f$ is constant through the time interval, then the analytic solution is:
   * \f[ \Psi(t) = \Psi_0 \exp\left( -\frac{V}{L} t \right) + \Psi_{ss} \left( 1 - \exp\left(
   * - \frac{V}{L} t\right) \right).\f]
   * @param stateVarReference \f$ \Psi_0 \f$
   * @param timeIncrement \f$ t \f$
   * @param localSlipRate \f$ V \f$
   * @return \f$ \Psi(t) \f$
   */
  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S stateVariableAt(const FrictionLawContext<Cfg, Type>& ctx,
                                              std::uint32_t pointIndex,
                                              real stateVarReference,
                                              real timeIncrement,
                                              S localSlipRate) {
    using std::expm1;
    using std::fmax;
    using std::log;
    using std::pow;
    const auto ltsFace = ctx.ltsFace;
    const real localMuW = ctx.data->muW[ltsFace][pointIndex];
    const real localSrW = ctx.data->srW[ltsFace][pointIndex];
    const real localA = ctx.data->a[ltsFace][pointIndex];
    const real localSl0 = ctx.data->sl0[ltsFace][pointIndex];

    // low-velocity steady state friction coefficient
    const S lowVelocityFriction =
        fmax(S(static_cast<real>(0)),
             S(ctx.data->f0[ltsFace][pointIndex]) -
                 S(ctx.data->b[ltsFace][pointIndex] - localA) *
                     log(localSlipRate / S(ctx.data->drParameters.rsSr0)));
    const S steadyStateFrictionCoefficient =
        S(localMuW) +
        (lowVelocityFriction - S(localMuW)) /
            pow(S(static_cast<real>(1.0)) + misc::power<8>(localSlipRate / S(localSrW)),
                static_cast<real>(1.0 / 8.0));
    // TODO: check again, if double precision is necessary here (earlier, there were cancellation
    // issues)
    // stated over V / (2 V_0) rather than under 2 V_0 / V: the reciprocal's derivative is
    // -2 V_0 / V^2, which at the slip-rate floor is 2e64 and leaves single precision, while the
    // state variable it belongs to is perfectly well behaved there
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
    const S localStateVariable =
        S(stateVarReference) + (steadyStateStateVariable - S(stateVarReference)) * exp1m;
    if constexpr (Type == common::RangeType::CPU) {
      assert((std::isfinite(valueOf(localStateVariable)) ||
              pointIndex >= misc::NumBoundaryGaussPoints<Cfg> * Cfg::NumSimulations) &&
             "Inf/NaN detected");
    }
    return localStateVariable;
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
    const auto stateVariable = stateVariableAt(
        ctx, pointIndex, stateVarReference, timeIncrement, dualCast<StateScalar>(slipRate));

    const S localStateVariable = dualCast<real>(stateVariable);
    const S localA = S(ctx.data->a[ctx.ltsFace][pointIndex]);
    const S cExpLog = localStateVariable / localA;
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
      const real localA = ctx.data->a[ctx.ltsFace][pointIndex];

      const real cLin = static_cast<real>(0.5) / ctx.data->drParameters.rsSr0;
      const real cExpLog = localStateVariable / localA;
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
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp
   * \left(\frac{\Psi}{a}\right)\right). \f]
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
   * Resample the state variable: its change in the step, so that the state variable lies in the
   * same polynomial space as the degrees of freedom (see resample).
   */
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void resampleStateVar(FrictionLawContext<Cfg, Type>& ctx) {
    const auto ltsFace = ctx.ltsFace;
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> deltaStateVar{};
    alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> resampledDeltaStateVar{};
    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      deltaStateVar[pointIndex] =
          ctx.stateVariableBuffer[pointIndex] - ctx.data->stateVariable[ltsFace][pointIndex];
    });

    resample(ctx, deltaStateVar, resampledDeltaStateVar);

    forEachPoint(ctx, [&](std::uint32_t pointIndex) {
      ctx.data->stateVariable[ltsFace][pointIndex] =
          ctx.data->stateVariable[ltsFace][pointIndex] + resampledDeltaStateVar[pointIndex];
    });
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FASTVELOCITYWEAKENINGLAW_H_
