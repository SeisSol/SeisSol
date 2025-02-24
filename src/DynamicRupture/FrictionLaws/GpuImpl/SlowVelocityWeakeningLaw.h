// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

namespace seissol::dr::friction_law::gpu {
template <class Derived, class TPMethod>
class SlowVelocityWeakeningLaw
    : public RateAndStateBase<SlowVelocityWeakeningLaw<Derived, TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<SlowVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {}

  // Note that we need double precision here, since single precision led to NaNs.
  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& ctx, double timeIncrement) {
    Derived::updateStateVariable(ctx, timeIncrement);
  }

  struct MuDetails {
    double a{};
    double c{};
    double ac{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext& ctx, double localStateVariable) {
    const double localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double log1 = std::log(ctx.data->drParameters.rsSr0 * localStateVariable / localSl0);
    const double c =
        0.5 / ctx.data->drParameters.rsSr0 *
        std::exp((ctx.data->drParameters.rsF0 + ctx.data->drParameters.rsB * log1) / localA);
    return MuDetails{localA, c, localA * c};
  }

  SEISSOL_DEVICE static double
      updateMu(FrictionLawContext& ctx, double localSlipRateMagnitude, MuDetails& details) {
    const double x = localSlipRateMagnitude * details.c;
    return details.a * std::asinh(x);
  }

  SEISSOL_DEVICE static double updateMuDerivative(FrictionLawContext& ctx,
                                                  double localSlipRateMagnitude,
                                                  MuDetails& details) {
    const double x = localSlipRateMagnitude * details.c;
    return details.ac / std::sqrt(std::pow(x, 2) + 1.0);
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
