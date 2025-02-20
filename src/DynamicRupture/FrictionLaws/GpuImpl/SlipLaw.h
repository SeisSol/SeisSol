// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::gpu {
template <typename TPMethod>
class SlipLaw : public SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::copyLtsTreeToLocal;

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {}

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
        localSl0 / localSlipRate * std::pow(localSlipRate * stateVarReference / localSl0, exp1);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
