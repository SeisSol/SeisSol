#ifndef SEISSOL_GPU_FRICTION_SOLVER_H
#define SEISSOL_GPU_FRICTION_SOLVER_H

#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include <omp.h>
#include <algorithm>

namespace seissol::dr::friction_law::gpu {

template <typename Derived>
class GpuFrictionSolver : public GpuBaseFrictionLaw {
  public:
  GpuFrictionSolver<Derived>(dr::DRParameters& drParameters) : GpuBaseFrictionLaw(drParameters) {}

  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture* dynRup,
                real fullUpdateTime,
                double timeWeights[CONVERGENCE_ORDER]) override {

    FrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);

    this->currLayerSize = layerData.getNumberOfCells();
    auto layerSize{this->currLayerSize};

    auto* impAndEta{this->impAndEta};
    auto* qInterpolatedPlus{this->qInterpolatedPlus};
    auto* qInterpolatedMinus{this->qInterpolatedMinus};
    auto* faultStresses{this->faultStresses};
    #pragma omp target teams loop          \
    map(from: faultStresses [0:layerSize]) \
    is_device_ptr(qInterpolatedPlus,       \
                  qInterpolatedMinus,      \
                  impAndEta)               \
    device(diviceId)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      Common::precomputeStressFromQInterpolated(faultStresses[ltsFace],
                                                impAndEta[ltsFace],
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace]);
    }

    static_cast<Derived*>(this)->preHook(stateVariableBuffer);

    // loop over all dynamic rupture faces, in this LTS layer
    #pragma omp parallel for schedule(static)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {

      // loop over sub time steps (i.e. quadrature points in time)
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
        static_cast<Derived*>(this)->updateFrictionAndSlip(faultStresses[ltsFace],
                                                           tractionResults[ltsFace],
                                                           stateVariableBuffer,
                                                           strengthBuffer,
                                                           ltsFace,
                                                           timeIndex);
      }
    }

    static_cast<Derived*>(this)->postHook(stateVariableBuffer);

    auto* ruptureTimePending{this->ruptureTimePending};
    auto* slipRateMagnitude{this->slipRateMagnitude};
    auto* ruptureTime{this->ruptureTime};
    #pragma omp target teams loop     \
    is_device_ptr(ruptureTimePending, \
                  slipRateMagnitude,  \
                  ruptureTime)        \
    firstprivate(mFullUpdateTime)     \
    device(diviceId)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      // output rupture front
      Common::saveRuptureFrontOutput(ruptureTimePending[ltsFace],
                                     ruptureTime[ltsFace],
                                     slipRateMagnitude[ltsFace],
                                     fullUpdateTime);
    }

    static_cast<Derived*>(this)->saveDynamicStressOutput();

    auto* peakSlipRate{this->peakSlipRate};
    auto* imposedStatePlus{this->imposedStatePlus};
    auto* imposedStateMinus{this->imposedStateMinus};
    auto* tractionResults{this->tractionResults};
    #pragma omp target teams loop             \
    map(to: faultStresses [0:layerSize],      \
            tractionResults[0:layerSize],     \
            timeWeights[0:CONVERGENCE_ORDER]) \
    is_device_ptr(peakSlipRate,       \
                  slipRateMagnitude,  \
                  imposedStatePlus,   \
                  imposedStateMinus,  \
                  qInterpolatedPlus,  \
                  qInterpolatedMinus, \
                  impAndEta)          \
    device(diviceId)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      Common::savePeakSlipRateOutput(slipRateMagnitude[ltsFace], peakSlipRate[ltsFace]);
      Common::postcomputeImposedStateFromNewStress(faultStresses[ltsFace],
                                                   tractionResults[ltsFace],
                                                   impAndEta[ltsFace],
                                                   imposedStatePlus[ltsFace],
                                                   imposedStateMinus[ltsFace],
                                                   qInterpolatedPlus[ltsFace],
                                                   qInterpolatedMinus[ltsFace],
                                                   timeWeights);
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_FRICTION_SOLVER_H
