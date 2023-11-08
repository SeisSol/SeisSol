#ifndef SEISSOL_GPU_SLIPLAW_H
#define SEISSOL_GPU_SLIPLAW_H

#include "DynamicRupture/FrictionLaws/GpuImpl/SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::gpu {
template <typename TPMethod>
class SlipLaw : public SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::copyLtsTreeToLocal;

  #pragma omp declare target
  void updateStateVariable(double timeIncrement) {
    auto* devSl0{this->sl0};
    auto* devStateVarReference{this->initialVariables.stateVarReference};
    auto* devLocalSlipRate{this->initialVariables.localSlipRate};
    auto* devStateVariableBuffer{this->stateVariableBuffer};
    auto layerSize{this->currLayerSize};

    // #pragma omp distribute
    #pragma omp target distribute map(in: devSl0[0:layerSize], devLocalSlipRate[0:layerSize], devStateVarReference[0:layerSize], out: devStateVariableBuffer[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {

        const double localSl0 = devSl0[ltsFace][pointIndex];
        const double localSlipRate = devLocalSlipRate[ltsFace][pointIndex];
        const double exp1 = std::exp(-localSlipRate * (timeIncrement / localSl0));

        const double stateVarReference = devStateVarReference[ltsFace][pointIndex];
        devStateVariableBuffer[ltsFace][pointIndex] =
            localSl0 / localSlipRate *
            std::pow(localSlipRate * stateVarReference / localSl0, exp1);
      }
    }
  }
  #pragma omp end declare target
};

} // namespace seissol::dr::friction_law::gpu
#endif // SEISSOL_GPU_SLIPLAW_H
