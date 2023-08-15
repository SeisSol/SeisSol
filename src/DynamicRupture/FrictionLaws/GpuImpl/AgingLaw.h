#ifndef SEISSOL_GPU_AGINGLAW_H
#define SEISSOL_GPU_AGINGLAW_H

#include "DynamicRupture/FrictionLaws/GpuImpl/SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::gpu {

template <typename Config, class TPMethod>
class AgingLaw : public SlowVelocityWeakeningLaw<Config, AgingLaw<Config, TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<Config, AgingLaw<Config, TPMethod>, TPMethod>::
      SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<Config, AgingLaw<Config, TPMethod>, TPMethod>::copyLtsTreeToLocal;

  void updateStateVariable(double timeIncrement) {
    auto* devSl0{this->sl0};
    auto* devStateVarReference{this->initialVariables.stateVarReference};
    auto* devLocalSlipRate{this->initialVariables.localSlipRate};
    auto* devStateVariableBuffer{this->stateVariableBuffer};

    sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints<Config>},
                       {misc::numPaddedPoints<Config>}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        const double localSl0 = devSl0[ltsFace][pointIndex];
        const double localSlipRate = devLocalSlipRate[ltsFace][pointIndex];
        const double exp1 = sycl::exp(-localSlipRate * (timeIncrement / localSl0));

        const double stateVarReference = devStateVarReference[ltsFace][pointIndex];
        devStateVariableBuffer[ltsFace][pointIndex] =
            static_cast<RealT>(stateVarReference * exp1 + localSl0 / localSlipRate * (1.0 - exp1));
      });
    });
  }
};

} // namespace seissol::dr::friction_law::gpu
#endif // SEISSOL_GPU_AGINGLAW_H
