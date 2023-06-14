#ifndef SEISSOL_BASE_FRICTION_SOLVER_H
#define SEISSOL_BASE_FRICTION_SOLVER_H

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "Numerical_aux/SyclFunctions.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include <algorithm>

namespace seissol::dr::friction_law::gpu {

template <typename Derived>
class BaseFrictionSolver : public FrictionSolverDetails {
  public:
  explicit BaseFrictionSolver<Derived>(dr::DRParameters* drParameters)
      : FrictionSolverDetails(drParameters) {}
  ~BaseFrictionSolver<Derived>() = default;

  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture const* const dynRup,
                real fullUpdateTime,
                const double timeWeights[CONVERGENCE_ORDER]) override {

    FrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->currLayerSize = layerData.getNumberOfCells();

    size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(double);
    auto timeWeightsCopy = this->queue.memcpy(devTimeWeights, &timeWeights[0], requiredNumBytes);

    {
      constexpr common::RangeType gpuRangeType{common::RangeType::GPU};

      auto* devImpAndEta{this->impAndEta};
      auto* devQInterpolatedPlus{this->qInterpolatedPlus};
      auto* devQInterpolatedMinus{this->qInterpolatedMinus};
      auto* devFaultStresses{this->faultStresses};

      sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
      this->queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          const auto ltsFace = item.get_group().get_group_id(0);
          const auto pointIndex = item.get_local_id(0);

          common::precomputeStressFromQInterpolated<gpuRangeType>(devFaultStresses[ltsFace],
                                                                  devImpAndEta[ltsFace],
                                                                  devQInterpolatedPlus[ltsFace],
                                                                  devQInterpolatedMinus[ltsFace],
                                                                  pointIndex);
        });
      });

      static_cast<Derived*>(this)->preHook(stateVariableBuffer);
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; ++timeIndex) {
        const real t0{this->drParameters->t0};
        const real dt = deltaT[timeIndex];
        auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
        const auto* devNucleationStressInFaultCS{this->nucleationStressInFaultCS};

        this->queue.submit([&](sycl::handler& cgh) {
          if (timeIndex == 0) {
            cgh.depends_on(timeWeightsCopy);
          }
          cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
            auto ltsFace = item.get_group().get_group_id(0);
            auto pointIndex = item.get_local_id(0);

            using StdMath = seissol::functions::SyclStdFunctions;
            common::adjustInitialStress<gpuRangeType, StdMath>(
                devInitialStressInFaultCS[ltsFace],
                devNucleationStressInFaultCS[ltsFace],
                fullUpdateTime,
                t0,
                dt,
                pointIndex);
          });
        });

        static_cast<Derived*>(this)->updateFrictionAndSlip(timeIndex);
      }
      static_cast<Derived*>(this)->postHook(stateVariableBuffer);

      auto* devRuptureTimePending{this->ruptureTimePending};
      auto* devSlipRateMagnitude{this->slipRateMagnitude};
      auto* devRuptureTime{this->ruptureTime};

      this->queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          auto ltsFace = item.get_group().get_group_id(0);
          auto pointIndex = item.get_local_id(0);
          common::saveRuptureFrontOutput<gpuRangeType>(devRuptureTimePending[ltsFace],
                                                       devRuptureTime[ltsFace],
                                                       devSlipRateMagnitude[ltsFace],
                                                       fullUpdateTime,
                                                       pointIndex);
        });
      });

      static_cast<Derived*>(this)->saveDynamicStressOutput();

      auto* devPeakSlipRate{this->peakSlipRate};
      auto* devImposedStatePlus{this->imposedStatePlus};
      auto* devImposedStateMinus{this->imposedStateMinus};
      auto* devTractionResults{this->tractionResults};
      auto* devTimeWeights{this->devTimeWeights};
      auto* devSpaceWeights{this->devSpaceWeights};
      auto* devEnergyData{this->energyData};
      auto* devGodunovData{this->godunovData};

      auto isFrictionEnergyRequired{this->drParameters->isEnergyFrictionRequired};
      this->queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          auto ltsFace = item.get_group().get_group_id(0);
          auto pointIndex = item.get_local_id(0);

          common::savePeakSlipRateOutput<gpuRangeType>(
              devSlipRateMagnitude[ltsFace], devPeakSlipRate[ltsFace], pointIndex);

          common::postcomputeImposedStateFromNewStress<gpuRangeType>(devFaultStresses[ltsFace],
                                                                     devTractionResults[ltsFace],
                                                                     devImpAndEta[ltsFace],
                                                                     devImposedStatePlus[ltsFace],
                                                                     devImposedStateMinus[ltsFace],
                                                                     devQInterpolatedPlus[ltsFace],
                                                                     devQInterpolatedMinus[ltsFace],
                                                                     devTimeWeights,
                                                                     pointIndex);

          if (isFrictionEnergyRequired) {
            common::computeFrictionEnergy<gpuRangeType>(devEnergyData[ltsFace],
                                                        devQInterpolatedPlus[ltsFace],
                                                        devQInterpolatedMinus[ltsFace],
                                                        devImpAndEta[ltsFace],
                                                        devTimeWeights,
                                                        devSpaceWeights,
                                                        devGodunovData[ltsFace],
                                                        pointIndex);
          }
        });
      });
      queue.wait_and_throw();
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_BASE_FRICTION_SOLVER_H
