// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::gpu {
template <typename Cfg, typename TPMethod>
class SlipLaw : public SlowVelocityWeakeningLaw<Cfg, SlipLaw<Cfg, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;
  using SlowVelocityWeakeningLaw<Cfg, SlipLaw<Cfg, TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<Cfg, SlipLaw<Cfg, TPMethod>, TPMethod>::copyStorageToLocal;

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {}

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext<Cfg>& ctx,
                                                 double timeIncrement) {
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double localSlipRate = ctx.initialVariables.localSlipRate;
    const double exp1v = std::exp(-localSlipRate * (timeIncrement / localSl0));

    const double stateVarReference = ctx.initialVariables.stateVarReference;
    ctx.stateVariableBuffer =
        localSl0 / localSlipRate * std::pow(localSlipRate * stateVarReference / localSl0, exp1v);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
