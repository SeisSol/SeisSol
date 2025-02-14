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
  static void updateStateVariable(FrictionLawContext& ctx, double timeIncrement) {
    Derived::updateStateVariable(ctx, timeIncrement);
  }

  struct Details {
    decltype(SlowVelocityWeakeningLaw::a) a;
    decltype(SlowVelocityWeakeningLaw::sl0) sl0;
    decltype(seissol::initializer::parameters::DRParameters::rsSr0) rsSr0;
    decltype(seissol::initializer::parameters::DRParameters::rsF0) rsF0;
    decltype(seissol::initializer::parameters::DRParameters::rsB) rsB;
  };

  struct MuDetails {
    double a{};
    double c{};
    double ac{};
  };

  static MuDetails getMuDetails(FrictionLawContext& ctx, double localStateVariable) {
    const double localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double log1 = sycl::log(ctx.data->drParameters.rsSr0 * localStateVariable / localSl0);
    const double c =
        0.5 / ctx.data->drParameters.rsSr0 *
        sycl::exp((ctx.data->drParameters.rsF0 + ctx.data->drParameters.rsB * log1) / localA);
    return MuDetails{localA, c, localA * c};
  }

  static double
      updateMu(FrictionLawContext& ctx, double localSlipRateMagnitude, MuDetails& details) {
    const double x = localSlipRateMagnitude * details.c;
    return details.a * sycl::asinh(x);
  }

  static double updateMuDerivative(FrictionLawContext& ctx,
                                   double localSlipRateMagnitude,
                                   MuDetails& details) {
    const double x = localSlipRateMagnitude * details.c;
    return details.ac / sycl::sqrt(sycl::pown(x, 2) + 1.0);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws,
   * we just copy the buffer into the member variable.
   */
  static void resampleStateVar(FrictionLawContext& ctx) {
    auto* stateVariable{ctx.data->stateVariable};

    stateVariable[ctx.ltsFace][ctx.pointIndex] = ctx.stateVariableBuffer;
  }

  static void executeIfNotConverged() {}
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
