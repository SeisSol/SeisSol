// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
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
  using SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod>::copyLtsTreeToLocal;

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& ctx, double timeIncrement) {
    auto* devSl0{ctx.data->sl0};
    auto& devStateVarReference{ctx.initialVariables.stateVarReference};
    auto& devLocalSlipRate{ctx.initialVariables.localSlipRate};
    auto& devStateVariableBuffer{ctx.stateVariableBuffer};

    const double localSl0 = devSl0[ctx.ltsFace][ctx.pointIndex];
    const double localSlipRate = devLocalSlipRate;
    const double exp1 = std::exp(-localSlipRate * (timeIncrement / localSl0));

    const double stateVarReference = devStateVarReference;
    devStateVariableBuffer =
        static_cast<real>(stateVarReference * exp1 + localSl0 / localSlipRate * (1.0 - exp1));
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
