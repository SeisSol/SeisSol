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
  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& ctx, double timeIncrement) {
    Derived::updateStateVariable(ctx, timeIncrement);
  }

  struct MuDetails {
    real a{};
    real cLin{};
    real cExpLog{};
    real cExp{};
    real acLin{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext& ctx, real localStateVariable) {
    const real localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const real log1 = std::log(ctx.data->drParameters.rsSr0 * localStateVariable / localSl0);
    const real localF0 = ctx.data->f0[ctx.ltsFace][ctx.pointIndex];
    const real localB = ctx.data->b[ctx.ltsFace][ctx.pointIndex];

    const real cLin = 0.5 / ctx.data->drParameters.rsSr0;
    const real cExpLog = (localF0 + localB * log1) / localA;
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

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws,
   * we just copy the buffer into the member variable.
   */
  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext& ctx) {
    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] = ctx.stateVariableBuffer;
  }

  SEISSOL_DEVICE static void executeIfNotConverged() {}
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
