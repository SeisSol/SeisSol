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

    Derived& self = *(static_cast<Derived*>(this));

    FrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);
    this->currLayerSize = layerData.getNumberOfCells();

    auto* queue{this->queue};

    size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(double);
    // auto timeWeightsCopy = this->queue.memcpy(devTimeWeights, &timeWeights[0], requiredNumBytes);
    double timeWeightsCopy[CONVERGENCE_ORDER];

    for (int i = 0; i < CONVERGENCE_ORDER; ++i) {
      timeWeightsCopy[i] = timeWeights[i];
      this->devTimeWeights[i] = timeWeights[i];
    }

     // map(to:timeWeightsCopy[0:CONVERGENCE_ORDER])

    // #pragma omp target teams
    {
      constexpr common::RangeType gpuRangeType{common::RangeType::GPU};

      auto* devImpAndEta{this->impAndEta};
      auto* devImpedanceMatrices{this->impedanceMatrices};
      auto* devQInterpolatedPlus{this->qInterpolatedPlus};
      auto* devQInterpolatedMinus{this->qInterpolatedMinus};
      auto* devFaultStresses{this->faultStresses};
      auto layerSize{this->currLayerSize};
      auto chunksize{this->chunksize};

      for (int chunk = 0; chunk < this->chunkcount; ++chunk)
      #pragma omp target depend(inout: queue[chunk]) device(TARGETDART_ANY) map(to:chunksize) map(to: CCHUNK(devImpAndEta), CCHUNK(devImpedanceMatrices), CCHUNK(devQInterpolatedPlus), CCHUNK(devQInterpolatedMinus)) map(from: CCHUNK(devFaultStresses)) nowait
      #pragma omp metadirective when(device={kind(nohost)}: teams distribute) default(parallel for)
      for (int ltsFace = 0; ltsFace < chunksize; ++ltsFace) {
        #pragma omp metadirective when(device={kind(nohost)}: parallel for) default(simd)
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
      self.preHook(devStateVariableBuffer);
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; ++timeIndex) {
        const real t0{this->drParameters->t0};
        const real dt = deltaT[timeIndex];
        auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
        const auto* devNucleationStressInFaultCS{this->nucleationStressInFaultCS};
        auto* devInitialPressure{this->initialPressure};
        const auto* devNucleationPressure{this->nucleationPressure};

        for (int chunk = 0; chunk < this->chunkcount; ++chunk)
        #pragma omp target depend(inout: queue[chunk]) device(TARGETDART_ANY) map(to:chunksize, timeIndex, dt, t0, fullUpdateTime) map(tofrom: CCHUNK(devInitialStressInFaultCS), CCHUNK(devInitialPressure)) map(to: CCHUNK(devNucleationStressInFaultCS), CCHUNK(devNucleationPressure)) nowait
        #pragma omp metadirective when(device={kind(nohost)}: teams distribute) default(parallel for)
        for (int ltsFace = 0; ltsFace < chunksize; ++ltsFace) {
          #pragma omp metadirective when(device={kind(nohost)}: parallel for) default(simd)
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

        self.updateFrictionAndSlip(timeIndex);
      }
      self.postHook(devStateVariableBuffer);

      auto* devRuptureTimePending{this->ruptureTimePending};
      auto* devSlipRateMagnitude{this->slipRateMagnitude};
      auto* devRuptureTime{this->ruptureTime};

      for (int chunk = 0; chunk < this->chunkcount; ++chunk)
      #pragma omp target depend(inout: queue[chunk]) device(TARGETDART_ANY) map(to:chunksize, fullUpdateTime) map(tofrom: CCHUNK(devRuptureTimePending)) map(from: CCHUNK(devRuptureTime)) map(to: CCHUNK(devSlipRateMagnitude)) nowait
      #pragma omp metadirective when(device={kind(nohost)}: teams distribute) default(parallel for)
      for (int ltsFace = 0; ltsFace < chunksize; ++ltsFace) {
        #pragma omp metadirective when(device={kind(nohost)}: parallel for) default(simd)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
          common::saveRuptureFrontOutput<gpuRangeType>(devRuptureTimePending[ltsFace],
                                                       devRuptureTime[ltsFace],
                                                       devSlipRateMagnitude[ltsFace],
                                                       fullUpdateTime,
                                                       pointIndex);
        }
      }

      self.saveDynamicStressOutput();

      auto* devPeakSlipRate{this->peakSlipRate};
      auto* devImposedStatePlus{this->imposedStatePlus};
      auto* devImposedStateMinus{this->imposedStateMinus};
      auto* devTractionResults{this->tractionResults};
      auto* devTimeWeights{this->devTimeWeights};
      auto* devSpaceWeights{this->devSpaceWeights};
      auto* devEnergyData{this->energyData};
      auto* devGodunovData{this->godunovData};

      /*
      #pragma omp target teams distribute depend(inout: queue[chunk]) device(TARGETDART_ANY) map(to: devTimeWeights[0:CONVERGENCE_ORDER], CCHUNK(devGodunovData), CCHUNK(devSlipRateMagnitude), CCHUNK(devFaultStresses), CCHUNK(devTractionResults), CCHUNK(devImpAndEta), CCHUNK(devImpedanceMatrices), CCHUNK(devQInterpolatedPlus), devQInterpolatedMinus) map(tofrom: CCHUNK(devPeakSlipRate), CCHUNK(devImposedStatePlus), CCHUNK(devImposedStateMinus), CCHUNK(devEnergyData)) nowait
      #pragma omp metadirective when( device={arch(nvptx)}: teams distribute) default(parallel for)
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp metadirective when( device={arch(nvptx)}: parallel for ) default(simd)
      */

      auto isFrictionEnergyRequired{this->drParameters->isFrictionEnergyRequired};
      
      for (int chunk = 0; chunk < this->chunkcount; ++chunk)
      #pragma omp target depend(inout: queue[chunk]) device(TARGETDART_ANY) map(to: chunksize, pointIndex) map(to: devTimeWeights[0:CONVERGENCE_ORDER], devSpaceWeights[0:misc::numPaddedPoints], CCHUNK(devGodunovData), CCHUNK(devSlipRateMagnitude), CCHUNK(devFaultStresses), CCHUNK(devTractionResults), CCHUNK(devImpAndEta), CCHUNK(devImpedanceMatrices), CCHUNK(devQInterpolatedPlus), CCHUNK(devQInterpolatedMinus)) map(tofrom: CCHUNK(devPeakSlipRate), CCHUNK(devImposedStatePlus), CCHUNK(devImposedStateMinus), CCHUNK(devEnergyData)) nowait
      #pragma omp metadirective when(device={kind(nohost)}: teams distribute) default(parallel for)
      for (int ltsFace = 0; ltsFace < chunksize; ++ltsFace) {
        #pragma omp metadirective when(device={kind(nohost)}: parallel for) default(simd)
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

      #pragma omp taskwait
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_BASE_FRICTION_SOLVER_H
