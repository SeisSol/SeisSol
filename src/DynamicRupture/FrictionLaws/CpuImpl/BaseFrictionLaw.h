// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_

#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Monitoring/Instrumentation.h"

#include <yaml-cpp/yaml.h>

namespace seissol::dr::friction_law::cpu {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Cfg, typename Derived>
class BaseFrictionLaw : public FrictionSolverImpl<Cfg> {
  private:
  size_t currLayerSize_{};

  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  explicit BaseFrictionLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverImpl<Cfg>(drParameters) {}

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& /*runtime*/) override {
    this->currLayerSize_ = layerData.size();
    FrictionSolverImpl<Cfg>::copyStorageToLocal(layerData);
    static_cast<Derived*>(this)->copyStorageToLocal(layerData);
  }

  /**
   * evaluates the current friction model
   */
  void evaluate(double fullUpdateTime,
                const FrictionSolver::FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& /*runtime*/) override {
    if (this->currLayerSize_ == 0) {
      return;
    }

    if constexpr (model::MaterialOf<Cfg>::SupportsDR) {

      SCOREP_USER_REGION_DEFINE(myRegionHandle)
      std::copy_n(frictionTime.deltaT.begin(), frictionTime.deltaT.size(), this->deltaT_);
      this->fullUpdateTime_ = fullUpdateTime;
      std::array<real, misc::TimeSteps<Cfg>> localTimeWeights{};
      std::copy_n(timeWeights, localTimeWeights.size(), localTimeWeights.begin());

      // loop over all dynamic rupture faces, in this LTS layer
#pragma omp parallel for schedule(static)
      for (std::size_t ltsFace = 0; ltsFace < this->currLayerSize_; ++ltsFace) {
        alignas(Alignment) ImposedState<Cfg, Executor::Host> imposedState{};
        alignas(Alignment) FaultStresses<Cfg, Executor::Host> faultStresses{};
        alignas(Alignment) FaultStresses<Cfg, Executor::Host> initialStress{};
        const auto etaPDamp = this->drParameters_.etaDampEnd > this->fullUpdateTime_
                                  ? this->drParameters_.etaDamp
                                  : static_cast<real>(1.0);

        SCOREP_USER_REGION_BEGIN(
            myRegionHandle, "computeDynamicRupturePreHook", SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRupturePreHook");
        // define some temporary variables
        std::array<real, misc::NumPaddedPoints<Cfg>> stateVariableBuffer{};
        std::array<real, misc::NumPaddedPoints<Cfg>> strengthBuffer{};

        static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
        LIKWID_MARKER_STOP("computeDynamicRupturePreHook");
        SCOREP_USER_REGION_END(myRegionHandle)

        // NOTE: this region now covers the whole sub time step pipeline -- the stress precompute
        // and the imposed state accumulation moved into the loop below and are no longer measured
        // on their own. Timings are therefore not comparable with the ones of the two regions
        // that used to surround it.
        SCOREP_USER_REGION_BEGIN(
            myRegionHandle, "computeDynamicRuptureTimeStepLoop", SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRuptureTimeStepLoop");
        TractionResults<Cfg, Executor::Host> tractionResults{};

        // loop over sub time steps (i.e. quadrature points in time
        real startTime = 0;
        real updateTime = this->fullUpdateTime_;
        for (std::size_t timeIndex = 0; timeIndex < misc::TimeSteps<Cfg>; timeIndex++) {
          startTime = updateTime;
          updateTime += this->deltaT_[timeIndex];

          common::precomputeStressFromQInterpolated<Cfg>(faultStresses,
                                                         this->impAndEta_[ltsFace],
                                                         this->impedanceMatrices_[ltsFace],
                                                         this->qInterpolatedPlus_[ltsFace],
                                                         this->qInterpolatedMinus_[ltsFace],
                                                         etaPDamp,
                                                         timeIndex);

          common::initializeTractionResults<Cfg>(faultStresses, tractionResults);

          const auto sourceCount = this->drParameters_.sourceCount;
          common::computeInitialStress<Cfg>(initialStress,
                                            &this->stressSourceInFaultCS_[ltsFace * sourceCount],
                                            &this->stressSourcePressure_[ltsFace * sourceCount],
                                            &this->stressSourceOnset_[ltsFace * sourceCount],
                                            &this->stressSourceRiseTime_[ltsFace * sourceCount],
                                            sourceCount,
                                            updateTime);

          static_cast<Derived*>(this)->updateFrictionAndSlip(faultStresses,
                                                             initialStress,
                                                             tractionResults,
                                                             stateVariableBuffer,
                                                             strengthBuffer,
                                                             ltsFace,
                                                             timeIndex);

          // time-dependent outputs
          common::saveRuptureFrontOutput<Cfg>(this->ruptureTimePending_[ltsFace],
                                              this->ruptureTime_[ltsFace],
                                              this->slipRateMagnitude_[ltsFace],
                                              startTime);

          static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace, startTime);

          common::savePeakSlipRateOutput<Cfg>(this->slipRateMagnitude_[ltsFace],
                                              this->peakSlipRate_[ltsFace]);

          if (this->drParameters_.isFrictionEnergyRequired &&
              this->drParameters_.isCheckAbortCriteraEnabled) {
            common::updateTimeSinceSlipRateBelowThreshold<Cfg>(
                this->slipRateMagnitude_[ltsFace],
                this->ruptureTimePending_[ltsFace],
                this->energyData_[ltsFace],
                this->deltaT_[timeIndex],
                this->drParameters_.terminatorSlipRateThreshold);
          }

          common::postcomputeImposedStateFromNewStress<Cfg>(imposedState,
                                                            faultStresses,
                                                            tractionResults,
                                                            this->impAndEta_[ltsFace],
                                                            this->impedanceMatrices_[ltsFace],
                                                            this->qInterpolatedPlus_[ltsFace],
                                                            this->qInterpolatedMinus_[ltsFace],
                                                            timeIndex,
                                                            localTimeWeights[timeIndex]);
        }
        LIKWID_MARKER_STOP("computeDynamicRuptureTimeStepLoop");
        SCOREP_USER_REGION_END(myRegionHandle)

        SCOREP_USER_REGION_BEGIN(
            myRegionHandle, "computeDynamicRupturePostHook", SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRupturePostHook");
        static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);

        LIKWID_MARKER_STOP("computeDynamicRupturePostHook");
        SCOREP_USER_REGION_END(myRegionHandle)

        SCOREP_USER_REGION_BEGIN(myRegionHandle,
                                 "computeDynamicRuptureFinalizeImposedState",
                                 SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRuptureFinalizeImposedState");
        common::finalizeImposedState<Cfg>(
            imposedState, this->imposedStatePlus_[ltsFace], this->imposedStateMinus_[ltsFace]);
        LIKWID_MARKER_STOP("computeDynamicRuptureFinalizeImposedState");
        SCOREP_USER_REGION_END(myRegionHandle)

        if (this->drParameters_.isFrictionEnergyRequired) {
          common::computeFrictionEnergy<Cfg>(this->energyData_[ltsFace],
                                             this->qInterpolatedPlus_[ltsFace],
                                             this->qInterpolatedMinus_[ltsFace],
                                             this->impAndEta_[ltsFace],
                                             localTimeWeights.data(),
                                             this->spaceWeights_,
                                             this->godunovData_[ltsFace],
                                             this->slipRateMagnitude_[ltsFace],
                                             this->drParameters_.energiesFromAcrossFaultVelocities);
        }
      }
    } else {
      logError() << "The material" << model::MaterialOf<Cfg>::Text
                 << "does not support DR friction law computations.";
    }
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
