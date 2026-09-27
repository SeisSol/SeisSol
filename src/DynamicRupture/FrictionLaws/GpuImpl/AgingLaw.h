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
    const S exp1v = exp(preexp1);

    const real stateVarReference = ctx.initialVariables.stateVarReference;
    // (L / V) (1 - exp(-V t / L)) is t times the mean of the relaxation over the step, which keeps
    // L / V and its derivative -L / V^2 out of the expression
    return S(stateVarReference) * exp1v + S(timeIncrement) * rs::relaxationWeight(-preexp1);
  }

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& __restrict ctx,
                                                 real timeIncrement) {
    ctx.stateVariableBuffer =
        stateVariableAt<real>(ctx, ctx.initialVariables.localSlipRate, timeIncrement);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
