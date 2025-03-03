// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

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

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {}

  // Note that we need double precision here, since single precision led to NaNs.
  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& ctx, double timeIncrement) {
    const double muW{ctx.data->drParameters.muW};

    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double localSlipRate = ctx.initialVariables.localSlipRate;

    const double steadyStateStateVariable = localSlipRate * localSl0 / ctx.data->drParameters.rsSr0;

    const double preexp1 = -ctx.data->drParameters.rsSr0 * (timeIncrement / localSl0);
    const double exp1 = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    const double localStateVariable =
        steadyStateStateVariable * exp1m + exp1 * ctx.initialVariables.stateVarReference;

    ctx.stateVariableBuffer = localStateVariable;
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

  struct MuDetails {
    double a{};
    double c{};
  };

  SEISSOL_DEVICE static MuDetails getMuDetails(FrictionLawContext& ctx, double localStateVariable) {
    const double localA = ctx.data->a[ctx.ltsFace][ctx.pointIndex];
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double c =
        ctx.data->drParameters.rsB * localStateVariable / (localStateVariable + localSl0);
    return MuDetails{localA, c};
  }

  SEISSOL_DEVICE static double
      updateMu(FrictionLawContext& ctx, double localSlipRateMagnitude, MuDetails& details) {
    return ctx.data->drParameters.rsF0 +
           details.a * localSlipRateMagnitude /
               (localSlipRateMagnitude + ctx.data->drParameters.rsSr0) -
           details.c;
  }

  SEISSOL_DEVICE static double updateMuDerivative(FrictionLawContext& ctx,
                                                  double localSlipRateMagnitude,
                                                  MuDetails& details) {
    // note that: d/dx (x/(x+c)) = ((x+c)-x)/(x+c)**2 = c/(x+c)**2
    const double divisor = (localSlipRateMagnitude + ctx.data->drParameters.rsSr0);
    return details.a * ctx.data->drParameters.rsSr0 / (divisor * divisor);
  }

  // no resampling
  SEISSOL_DEVICE static void resampleStateVar(FrictionLawContext& ctx) {
    ctx.data->stateVariable[ctx.ltsFace][ctx.pointIndex] = ctx.stateVariableBuffer;
  }

  SEISSOL_DEVICE static void executeIfNotConverged() {}
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
