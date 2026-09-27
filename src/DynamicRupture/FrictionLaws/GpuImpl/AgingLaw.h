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

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& __restrict ctx,
                                                 double timeIncrement) {
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const real localSlipRate = ctx.initialVariables.localSlipRate;
    const double preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const double exp1v = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);

    // the relaxation towards L / V with a single exponential, which keeps L / V as its exact fixed
    // point (see FastVelocityWeakeningLaw::updateStateVariable); with exp once the step relaxes
    // more than half the way, since mu takes the logarithm of the state and the form with expm1
    // alone may then round it to zero
    const double stateVarReference = ctx.initialVariables.stateVarReference;
    const double steadyStateStateVariable = static_cast<double>(localSl0) / localSlipRate;
    ctx.stateVariableBuffer = static_cast<real>(
        exp1m < 0.5
            ? stateVarReference + (steadyStateStateVariable - stateVarReference) * exp1m
            : steadyStateStateVariable + (stateVarReference - steadyStateStateVariable) * exp1v);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
