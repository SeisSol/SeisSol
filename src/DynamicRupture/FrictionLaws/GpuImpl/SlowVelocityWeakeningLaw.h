// SPDX-FileCopyrightText: 2022 SeisSol Group
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
template <typename Cfg, class Derived, class TPMethod>
class SlowVelocityWeakeningLaw
    : public RateAndStateBase<Cfg, SlowVelocityWeakeningLaw<Cfg, Derived, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;
  using RateAndStateBase<Cfg, SlowVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  // Note that we need double precision here, since single precision led to NaNs.
  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext<Cfg>& ctx,
                                                 double timeIncrement) {
    Derived::updateStateVariable(ctx, timeIncrement);
  }

  struct MuDetails {
    double a{};
    double c{};
    double ac{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext<Cfg>& ctx,
                                               double localStateVariable) {
    const double localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double log1 = std::log(ctx.data->drParameters.rsSr0 * localStateVariable / localSl0);
    const double c =
        0.5 / ctx.data->drParameters.rsSr0 *
        std::exp((ctx.data->drParameters.rsF0 + ctx.data->drParameters.rsB * log1) / localA);
    return MuDetails{localA, c, localA * c};
  }

  SEISSOL_DEVICE static double updateMu(FrictionLawContext<Cfg>& ctx,
                                        double localSlipRateMagnitude,
                                        const MuDetails& details) {
    const double x = localSlipRateMagnitude * details.c;
    return details.a * std::asinh(x);
  }

  SEISSOL_DEVICE static double updateMuDerivative(FrictionLawContext<Cfg>& ctx,
                                                  double localSlipRateMagnitude,
                                                  const MuDetails& details) {
    const double x = localSlipRateMagnitude * details.c;
    return details.ac / std::sqrt(x * x + 1.0);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws,
   * we just copy the buffer into the member variable.
   */
  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext<Cfg>& ctx) {
    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] = ctx.stateVariableBuffer;
  }

  SEISSOL_DEVICE static void executeIfNotConverged() {}
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
