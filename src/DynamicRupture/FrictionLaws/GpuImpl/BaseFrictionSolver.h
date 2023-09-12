#ifndef SEISSOL_BASE_FRICTION_SOLVER_H
#define SEISSOL_BASE_FRICTION_SOLVER_H

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
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
    // auto timeWeightsCopy = this->queue.memcpy(devTimeWeights, &timeWeights[0], requiredNumBytes);
    double timeWeightsCopy[CONVERGENCE_ORDER];

    for (int i = 0; i < CONVERGENCE_ORDER; ++i) {
      timeWeightsCopy[i] = timeWeights[i];
      this->devTimeWeights[i] = timeWeights[i];
    }

    #pragma omp target teams // map(to:timeWeightsCopy[0:CONVERGENCE_ORDER])
    {
      constexpr common::RangeType gpuRangeType{common::RangeType::GPU};

      auto* devImpAndEta{this->impAndEta};
      auto* devImpedanceMatrices{this->impedanceMatrices};
      auto* devQInterpolatedPlus{this->qInterpolatedPlus};
      auto* devQInterpolatedMinus{this->qInterpolatedMinus};
      auto* devFaultStresses{this->faultStresses};

      #pragma omp distribute
      for (int ltsFace = 0; ltsFace < this->currLayerSize; ++ltsFace) {
        #pragma omp parallel for
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
          common::precomputeStressFromQInterpolated<gpuRangeType>(devFaultStresses[ltsFace],
                                                                  devImpAndEta[ltsFace],
                                                                  devImpedanceMatrices[ltsFace],
                                                                  devQInterpolatedPlus[ltsFace],
                                                                  devQInterpolatedMinus[ltsFace],
                                                                  pointIndex);
        }
      }

      auto* devStateVariableBuffer = this->stateVariableBuffer;
      static_cast<Derived*>(this)->preHook(devStateVariableBuffer);
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; ++timeIndex) {
        const real t0{this->drParameters->t0};
        const real dt = deltaT[timeIndex];
        auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
        const auto* devNucleationStressInFaultCS{this->nucleationStressInFaultCS};
        auto* devInitialPressure{this->initialPressure};
        const auto* devNucleationPressure{this->nucleationPressure};

        #pragma omp distribute
        for (int ltsFace = 0; ltsFace < this->currLayerSize; ++ltsFace) {
          #pragma omp parallel for
          for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
          // if (timeIndex == 0) {cgh.depends_on(timeWeightsCopy);}
            common::adjustInitialStress<gpuRangeType>(
                devInitialStressInFaultCS[ltsFace],
                devNucleationStressInFaultCS[ltsFace],
                devInitialPressure[ltsFace],
                devNucleationPressure[ltsFace],
                fullUpdateTime,
                t0,
                dt,
                pointIndex);
          }
        }

        static_cast<Derived*>(this)->updateFrictionAndSlip(timeIndex);
      }
      static_cast<Derived*>(this)->postHook(devStateVariableBuffer);

      auto* devRuptureTimePending{this->ruptureTimePending};
      auto* devSlipRateMagnitude{this->slipRateMagnitude};
      auto* devRuptureTime{this->ruptureTime};

      #pragma omp distribute
      for (int ltsFace = 0; ltsFace < this->currLayerSize; ++ltsFace) {
        #pragma omp parallel for
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
          common::saveRuptureFrontOutput<gpuRangeType>(devRuptureTimePending[ltsFace],
                                                       devRuptureTime[ltsFace],
                                                       devSlipRateMagnitude[ltsFace],
                                                       fullUpdateTime,
                                                       pointIndex);
        }
      }

      static_cast<Derived*>(this)->saveDynamicStressOutput();

      auto* devPeakSlipRate{this->peakSlipRate};
      auto* devImposedStatePlus{this->imposedStatePlus};
      auto* devImposedStateMinus{this->imposedStateMinus};
      auto* devTractionResults{this->tractionResults};
      auto* devTimeWeights{this->devTimeWeights};
      auto* devSpaceWeights{this->devSpaceWeights};
      auto* devEnergyData{this->energyData};
      auto* devGodunovData{this->godunovData};

      auto isFrictionEnergyRequired{this->drParameters->isFrictionEnergyRequired};
      #pragma omp distribute
      for (int ltsFace = 0; ltsFace < this->currLayerSize; ++ltsFace) {
        #pragma omp parallel for
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {

          common::savePeakSlipRateOutput<gpuRangeType>(
              devSlipRateMagnitude[ltsFace], devPeakSlipRate[ltsFace], pointIndex);

          common::postcomputeImposedStateFromNewStress<gpuRangeType>(devFaultStresses[ltsFace],
                                                                     devTractionResults[ltsFace],
                                                                     devImpAndEta[ltsFace],
                                                                     devImpedanceMatrices[ltsFace],
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
        }
      }
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_BASE_FRICTION_SOLVER_H
