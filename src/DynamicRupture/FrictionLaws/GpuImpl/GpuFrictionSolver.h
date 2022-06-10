#ifndef SEISSOL_GPU_FRICTION_SOLVER_H
#define SEISSOL_GPU_FRICTION_SOLVER_H

#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "Numerical_aux/SyclFunctions.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include <algorithm>

namespace seissol::dr::friction_law::gpu {

template <typename Derived>
class GpuFrictionSolver : public GpuBaseFrictionLaw {
  public:
  GpuFrictionSolver<Derived>(dr::DRParameters* drParameters) : GpuBaseFrictionLaw(drParameters) {}

  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture const* const dynRup,
                real fullUpdateTime,
                const double timeWeights[CONVERGENCE_ORDER]) override {

    FrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->currLayerSize = layerData.getNumberOfCells();

    size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(double);
    this->queue.memcpy(devTimeWeights, &timeWeights[0], requiredNumBytes).wait();

    {
      constexpr common::RangeType gpuRangeType{common::RangeType::GPU};
      auto layerSize{this->currLayerSize};

      auto* impAndEta{this->impAndEta};
      auto* impedanceMatrices{this->impedanceMatrices};
      auto* qInterpolatedPlus{this->qInterpolatedPlus};
      auto* qInterpolatedMinus{this->qInterpolatedMinus};
      auto* faultStresses{this->faultStresses};

      sycl::nd_range rng{{layerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
      this->queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          const auto ltsFace = item.get_group().get_group_id(0);
          const auto pointIndex = item.get_local_id(0);

          common::precomputeStressFromQInterpolated<gpuRangeType>(faultStresses[ltsFace],
                                                                  impedanceMatrices[ltsFace],
                                                                  qInterpolatedPlus[ltsFace],
                                                                  qInterpolatedMinus[ltsFace],
                                                                  pointIndex);
        });
      });

      static_cast<Derived*>(this)->preHook(stateVariableBuffer);
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; ++timeIndex) {
        const real t0{this->drParameters->t0};
        const real dt = deltaT[timeIndex];
        auto* initialStressInFaultCS{this->initialStressInFaultCS};
        const auto* nucleationStressInFaultCS{this->nucleationStressInFaultCS};

        this->queue.submit([&](sycl::handler& cgh) {
          cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
            auto ltsFace = item.get_group().get_group_id(0);
            auto pointIndex = item.get_local_id(0);

            using StdMath = seissol::functions::SyclStdFunctions;
            common::adjustInitialStress<gpuRangeType, StdMath>(initialStressInFaultCS[ltsFace],
                                                               nucleationStressInFaultCS[ltsFace],
                                                               fullUpdateTime,
                                                               t0,
                                                               dt,
                                                               pointIndex);
          });
        });

        static_cast<Derived*>(this)->updateFrictionAndSlip(timeIndex);
      }
      static_cast<Derived*>(this)->postHook(stateVariableBuffer);

      auto* ruptureTimePending{this->ruptureTimePending};
      auto* slipRateMagnitude{this->slipRateMagnitude};
      auto* ruptureTime{this->ruptureTime};

      this->queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          auto ltsFace = item.get_group().get_group_id(0);
          auto pointIndex = item.get_local_id(0);
          common::saveRuptureFrontOutput<gpuRangeType>(ruptureTimePending[ltsFace],
                                                       ruptureTime[ltsFace],
                                                       slipRateMagnitude[ltsFace],
                                                       fullUpdateTime,
                                                       pointIndex);
        });
      });

      static_cast<Derived*>(this)->saveDynamicStressOutput();

      auto* peakSlipRate{this->peakSlipRate};
      auto* imposedStatePlus{this->imposedStatePlus};
      auto* imposedStateMinus{this->imposedStateMinus};
      auto* tractionResults{this->tractionResults};
      auto* devTimeWeights{this->devTimeWeights};
      auto* devSpaceWeights{this->devSpaceWeights};
      auto* energyData{this->energyData};
      auto* godunovData{this->godunovData};

      this->queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          auto ltsFace = item.get_group().get_group_id(0);
          auto pointIndex = item.get_local_id(0);

          common::savePeakSlipRateOutput<gpuRangeType>(
              slipRateMagnitude[ltsFace], peakSlipRate[ltsFace], pointIndex);

          common::postcomputeImposedStateFromNewStress<gpuRangeType>(faultStresses[ltsFace],
                                                                     tractionResults[ltsFace],
                                                                     impedanceMatrices[ltsFace],
                                                                     imposedStatePlus[ltsFace],
                                                                     imposedStateMinus[ltsFace],
                                                                     qInterpolatedPlus[ltsFace],
                                                                     qInterpolatedMinus[ltsFace],
                                                                     devTimeWeights,
                                                                     pointIndex);

          common::computeFrictionEnergy<gpuRangeType>(energyData[ltsFace],
                                                      qInterpolatedPlus[ltsFace],
                                                      qInterpolatedMinus[ltsFace],
                                                      impedanceMatrices[ltsFace],
                                                      devTimeWeights,
                                                      devSpaceWeights,
                                                      godunovData[ltsFace],
                                                      pointIndex);
        });
      });
      queue.wait_and_throw();
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_FRICTION_SOLVER_H
