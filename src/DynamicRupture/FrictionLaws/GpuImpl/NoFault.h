#ifndef SEISSOL_GPU_NOFAULT_H
#define SEISSOL_GPU_NOFAULT_H

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"

namespace seissol::dr::friction_law::gpu {

class NoFault : public BaseFrictionSolver<NoFault> {
  public:
  NoFault(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionSolver<NoFault>(drParameters) {};

  void copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                      const seissol::initializer::DynamicRupture* const dynRup,
                                      real fullUpdateTime) override {}

  void updateFrictionAndSlip(unsigned timeIndex) {
    auto* devTraction1{this->traction1};
    auto* devTraction2{this->traction2};
    auto* devFaultStresses{this->faultStresses};
    auto* devTractionResults{this->tractionResults};

    auto chunksize{this->chunksize};
    auto layerSize{this->layerSize};
    auto* queue{this->queue};

    for (int chunk = 0; chunk < this->chunkcount; ++chunk)
    #pragma omp target depend(inout: queue[chunk]) device(TARGETDART_ANY) map(to:CCHUNK(devFaultStresses))) map(from: CCHUNK(devTractionResults), CCHUNK(traction1), CCHUNK(traction2)) nowait
    #pragma omp metadirective when(device={kind(nohost)}: teams distribute) default(parallel for)
      CCHUNKLOOP(ltsFace) {
        #pragma omp metadirective when(device={kind(nohost)}: parallel for) default(simd)
        for (int pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        auto& faultStresses = devFaultStresses[ltsFace];
        auto& tractionResults = devTractionResults[ltsFace];

        // calculate traction
        tractionResults.traction1[timeIndex][pointIndex] =
            faultStresses.traction1[timeIndex][pointIndex];
        tractionResults.traction2[timeIndex][pointIndex] =
            faultStresses.traction2[timeIndex][pointIndex];
        devTraction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
        devTraction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
      }
    }
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  void saveDynamicStressOutput() {}

  void preHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {}
  void postHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {}

  protected:
};

} // namespace seissol::dr::friction_law::gpu

#endif
