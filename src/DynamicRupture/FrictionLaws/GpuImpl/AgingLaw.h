// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::gpu {

template <class TPMethod>
class AgingLaw : public SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod>::copyStorageToLocal;

  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S>
  SEISSOL_DEVICE static S
      stateVariableAt(FrictionLawContext& __restrict ctx, S localSlipRate, real timeIncrement) {
    using std::exp;
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const S preexp1 = -localSlipRate * S(timeIncrement / localSl0);
    // (L / V) (1 - exp(-V t / L)) is t times the mean of the relaxation over the step, which keeps
    // L / V and its derivative -L / V^2 out of the expression
    const S weight = rs::relaxationWeight(-preexp1);

    const real stateVarReference = ctx.initialVariables.stateVarReference;
    // the relaxation towards L / V with a single exponential, psi0 + (L / V - psi0) (1 - exp(p)),
    // which keeps L / V as its fixed point (see FastVelocityWeakeningLaw::updateStateVariable);
    // with 1 - exp(p) = -p weight, it reads psi0 + weight (t + psi0 p) and leaves L / V out as
    // well. With exp once the step relaxes more than half the way, since mu takes the logarithm
    // of the state and the form with expm1 alone may then round it to zero.
    if (valueOf(-preexp1 * weight) < static_cast<real>(0.5)) {
      return S(stateVarReference) + weight * (S(timeIncrement) + S(stateVarReference) * preexp1);
    }
    return S(stateVarReference) * exp(preexp1) + S(timeIncrement) * weight;
  }

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& __restrict ctx,
                                                 real timeIncrement) {
    ctx.stateVariableBuffer =
        stateVariableAt<real>(ctx, ctx.initialVariables.localSlipRate, timeIncrement);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
