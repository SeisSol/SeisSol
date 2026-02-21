// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_

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
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  private:
  size_t currLayerSize_{};

  public:
  explicit BaseFrictionLaw(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolver(drParameters) {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& /*runtime*/) override {
    this->currLayerSize_ = layerData.size();
    BaseFrictionLaw::copyStorageToLocal(layerData);
    static_cast<Derived*>(this)->copyStorageToLocal(layerData);
  }

  /**
   * evaluates the current friction model
   */
  void evaluate(real fullUpdateTime,
                const FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& /*runtime*/) override {
    if (this->currLayerSize_ == 0) {
      return;
    }

    if constexpr (model::MaterialT::SupportsDR) {

      SCOREP_USER_REGION_DEFINE(myRegionHandle)
      std::copy_n(frictionTime.deltaT.begin(), frictionTime.deltaT.size(), this->deltaT_);
      this->mFullUpdateTime_ = fullUpdateTime;

      // loop over all dynamic rupture faces, in this LTS layer
#pragma omp parallel for schedule(static)
      for (std::size_t ltsFace = 0; ltsFace < this->currLayerSize_; ++ltsFace) {
        alignas(Alignment) FaultStresses<Executor::Host> faultStresses{};
        SCOREP_USER_REGION_BEGIN(
            myRegionHandle, "computeDynamicRupturePrecomputeStress", SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRupturePrecomputeStress");
        const auto etaPDamp =
            drParameters_->etaDampEnd > this->mFullUpdateTime_ ? drParameters_->etaDamp : 1.0;
        common::precomputeStressFromQInterpolated(faultStresses,
                                                  impAndEta_[ltsFace],
                                                  impedanceMatrices_[ltsFace],
                                                  qInterpolatedPlus_[ltsFace],
                                                  qInterpolatedMinus_[ltsFace],
                                                  etaPDamp);
        LIKWID_MARKER_STOP("computeDynamicRupturePrecomputeStress");
        SCOREP_USER_REGION_END(myRegionHandle)

        SCOREP_USER_REGION_BEGIN(
            myRegionHandle, "computeDynamicRupturePreHook", SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRupturePreHook");
        // define some temporary variables
        std::array<real, misc::NumPaddedPoints> stateVariableBuffer{0};
        std::array<real, misc::NumPaddedPoints> strengthBuffer{0};

        static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
        LIKWID_MARKER_STOP("computeDynamicRupturePreHook");
        SCOREP_USER_REGION_END(myRegionHandle)

        SCOREP_USER_REGION_BEGIN(myRegionHandle,
                                 "computeDynamicRuptureUpdateFrictionAndSlip",
                                 SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRuptureUpdateFrictionAndSlip");
        TractionResults<Executor::Host> tractionResults = {};

        // loop over sub time steps (i.e. quadrature points in time
        real startTime = 0;
        real updateTime = this->mFullUpdateTime_;
        for (std::size_t timeIndex = 0; timeIndex < misc::TimeSteps; timeIndex++) {
          startTime = updateTime;
          updateTime += this->deltaT_[timeIndex];
          for (unsigned i = 0; i < this->drParameters_->nucleationCount; ++i) {
            common::adjustInitialStress(
                initialStressInFaultCS_[ltsFace],
                nucleationStressInFaultCS_[ltsFace * this->drParameters_->nucleationCount + i],
                initialPressure_[ltsFace],
                nucleationPressure_[ltsFace * this->drParameters_->nucleationCount + i],
                updateTime,
                this->drParameters_->t0[i],
                this->drParameters_->s0[i],
                this->deltaT_[timeIndex]);
          }

          static_cast<Derived*>(this)->updateFrictionAndSlip(faultStresses,
                                                             tractionResults,
                                                             stateVariableBuffer,
                                                             strengthBuffer,
                                                             ltsFace,
                                                             timeIndex);

          // time-dependent outputs
          common::saveRuptureFrontOutput(ruptureTimePending_[ltsFace],
                                         ruptureTime_[ltsFace],
                                         slipRateMagnitude_[ltsFace],
                                         startTime);

          static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace, startTime);

          common::savePeakSlipRateOutput(slipRateMagnitude_[ltsFace], peakSlipRate_[ltsFace]);

          if (this->drParameters_->isFrictionEnergyRequired &&
              this->drParameters_->isCheckAbortCriteraEnabled) {
            common::updateTimeSinceSlipRateBelowThreshold(
                slipRateMagnitude_[ltsFace],
                ruptureTimePending_[ltsFace],
                energyData_[ltsFace],
                this->deltaT_[timeIndex],
                this->drParameters_->terminatorSlipRateThreshold);
          }
        }
        LIKWID_MARKER_STOP("computeDynamicRuptureUpdateFrictionAndSlip");
        SCOREP_USER_REGION_END(myRegionHandle)

        SCOREP_USER_REGION_BEGIN(
            myRegionHandle, "computeDynamicRupturePostHook", SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRupturePostHook");
        static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);

        LIKWID_MARKER_STOP("computeDynamicRupturePostHook");
        SCOREP_USER_REGION_END(myRegionHandle)

        SCOREP_USER_REGION_BEGIN(myRegionHandle,
                                 "computeDynamicRupturePostcomputeImposedState",
                                 SCOREP_USER_REGION_TYPE_COMMON)
        LIKWID_MARKER_START("computeDynamicRupturePostcomputeImposedState");
        common::postcomputeImposedStateFromNewStress(faultStresses,
                                                     tractionResults,
                                                     impAndEta_[ltsFace],
                                                     impedanceMatrices_[ltsFace],
                                                     imposedStatePlus_[ltsFace],
                                                     imposedStateMinus_[ltsFace],
                                                     qInterpolatedPlus_[ltsFace],
                                                     qInterpolatedMinus_[ltsFace],
                                                     timeWeights);
        LIKWID_MARKER_STOP("computeDynamicRupturePostcomputeImposedState");
        SCOREP_USER_REGION_END(myRegionHandle)

        if (this->drParameters_->isFrictionEnergyRequired) {
          common::computeFrictionEnergy(energyData_[ltsFace],
                                        qInterpolatedPlus_[ltsFace],
                                        qInterpolatedMinus_[ltsFace],
                                        impAndEta_[ltsFace],
                                        timeWeights,
                                        spaceWeights_,
                                        godunovData_[ltsFace],
                                        slipRateMagnitude_[ltsFace],
                                        this->drParameters_->energiesFromAcrossFaultVelocities);
        }
      }
    } else {
      logError() << "The material" << model::MaterialT::Text
                 << "does not support DR friction law computations.";
    }
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
