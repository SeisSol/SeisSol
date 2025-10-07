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

template <typename Cfg, class TPMethod>
class AgingLaw : public SlowVelocityWeakeningLaw<Cfg, AgingLaw<Cfg, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;
  using SlowVelocityWeakeningLaw<Cfg, AgingLaw<Cfg, TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<Cfg, AgingLaw<Cfg, TPMethod>, TPMethod>::copyStorageToLocal;

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext<Cfg>& ctx,
                                                 double timeIncrement) {
    const double localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const double localSlipRate = ctx.initialVariables.localSlipRate;
    const double preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const double exp1v = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);

    const double stateVarReference = ctx.initialVariables.stateVarReference;
    ctx.stateVariableBuffer =
        static_cast<real>(stateVarReference * exp1v + localSl0 / localSlipRate * exp1m);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
