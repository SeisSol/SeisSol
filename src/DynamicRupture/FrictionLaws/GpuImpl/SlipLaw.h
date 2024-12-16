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
    auto* queue{this->queue};

    auto chunksize{this->chunksize};

    for (int chunk = 0; chunk < this->chunkcount; ++chunk)
    #pragma omp target depend(inout: queue[0]) device(TARGETDART_ANY) map(to:chunksize, timeIncrement) map(to: CCHUNK(devSl0), CCHUNK(devLocalSlipRate), CCHUNK(devStateVarReference)) map(from: CCHUNK(devStateVariableBuffer)) nowait
    #pragma omp metadirective when(device={kind(nohost)}: teams distribute) default(parallel for)
      CCHUNKLOOP(ltsFace) {
        #pragma omp metadirective when(device={kind(nohost)}: parallel for) default(simd)
        for (int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {

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
