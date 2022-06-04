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

    size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(double);
    omp_target_memcpy(devTimeWeights, &timeWeights[0], requiredNumBytes, 0, 0, deviceId, hostId);

    requiredNumBytes = CONVERGENCE_ORDER * sizeof(real);
    omp_target_memcpy(devDeltaT, &deltaT[0], requiredNumBytes, 0, 0, deviceId, hostId);

    {
      auto layerSize{this->currLayerSize};

      auto* impAndEta{this->impAndEta};
      auto* qInterpolatedPlus{this->qInterpolatedPlus};
      auto* qInterpolatedMinus{this->qInterpolatedMinus};
      auto* faultStresses{this->faultStresses};

      #pragma omp target teams loop \
      is_device_ptr(faultStresses, qInterpolatedPlus, qInterpolatedMinus, impAndEta) \
      device(deviceId) nowait
      for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        Common::precomputeStressFromQInterpolated(faultStresses[ltsFace],
                                                  impAndEta[ltsFace],
                                                  qInterpolatedPlus[ltsFace],
                                                  qInterpolatedMinus[ltsFace]);
      }

      static_cast<Derived*>(this)->preHook(stateVariableBuffer);
      static_cast<Derived*>(this)->updateFrictionAndSlip();
      static_cast<Derived*>(this)->postHook(stateVariableBuffer);

      auto* ruptureTimePending{this->ruptureTimePending};
      auto* slipRateMagnitude{this->slipRateMagnitude};
      auto* ruptureTime{this->ruptureTime};

      #pragma omp target teams loop \
      is_device_ptr(ruptureTimePending, slipRateMagnitude, ruptureTime) \
      device(deviceId) nowait
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
      auto* devTimeWeights{this->devTimeWeights};

      #pragma omp target teams loop     \
      is_device_ptr(faultStresses,      \
                    tractionResults,    \
                    peakSlipRate,       \
                    slipRateMagnitude,  \
                    imposedStatePlus,   \
                    imposedStateMinus,  \
                    qInterpolatedPlus,  \
                    qInterpolatedMinus, \
                    impAndEta,          \
                    devTimeWeights)     \
      device(deviceId) nowait
      for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        Common::savePeakSlipRateOutput(slipRateMagnitude[ltsFace], peakSlipRate[ltsFace]);
        Common::postcomputeImposedStateFromNewStress(faultStresses[ltsFace],
                                                     tractionResults[ltsFace],
                                                     impAndEta[ltsFace],
                                                     imposedStatePlus[ltsFace],
                                                     imposedStateMinus[ltsFace],
                                                     qInterpolatedPlus[ltsFace],
                                                     qInterpolatedMinus[ltsFace],
                                                     devTimeWeights);
      }

      #pragma omp taskwait
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_FRICTION_SOLVER_H
