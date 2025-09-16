// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_

#include <Memory/Descriptor/DynamicRupture.h>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Monitoring/Instrumentation.h"

namespace seissol::dr::friction_law::cpu {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  private:
  size_t currLayerSize{};

  public:
  explicit BaseFrictionLaw(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolver(drParameters) {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& runtime) override {
    this->currLayerSize = layerData.size();
    BaseFrictionLaw::copyStorageToLocal(layerData);
    static_cast<Derived*>(this)->copyStorageToLocal(layerData);
  }

  /**
   * evaluates the current friction model
   */
  void evaluate(real fullUpdateTime,
                const FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (this->currLayerSize == 0) {
      return;
    }

    SCOREP_USER_REGION_DEFINE(myRegionHandle)
    std::copy_n(frictionTime.deltaT.begin(), frictionTime.deltaT.size(), this->deltaT);
    this->sumDt = frictionTime.sumDt;
    this->mFullUpdateTime = fullUpdateTime;

    // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t ltsFace = 0; ltsFace < this->currLayerSize; ++ltsFace) {
      alignas(Alignment) FaultStresses<Executor::Host> faultStresses{};
      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePrecomputeStress", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePrecomputeStress");
      const auto etaPDamp =
          drParameters->etaStop > this->mFullUpdateTime ? drParameters->etaHack : 1.0;
      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                impedanceMatrices[ltsFace],
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace],
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
      real updateTime = this->mFullUpdateTime;
      for (std::size_t timeIndex = 0; timeIndex < misc::TimeSteps; timeIndex++) {
        startTime = updateTime;
        updateTime += this->deltaT[timeIndex];
        for (unsigned i = 0; i < this->drParameters->nucleationCount; ++i) {
          common::adjustInitialStress(
              initialStressInFaultCS[ltsFace],
              nucleationStressInFaultCS[ltsFace * this->drParameters->nucleationCount + i],
              initialPressure[ltsFace],
              nucleationPressure[ltsFace * this->drParameters->nucleationCount + i],
              updateTime,
              this->drParameters->t0[i],
              this->drParameters->s0[i],
              this->deltaT[timeIndex]);
        }

        static_cast<Derived*>(this)->updateFrictionAndSlip(faultStresses,
                                                           tractionResults,
                                                           stateVariableBuffer,
                                                           strengthBuffer,
                                                           ltsFace,
                                                           timeIndex);

        // time-dependent outputs
        common::saveRuptureFrontOutput(ruptureTimePending[ltsFace],
                                       ruptureTime[ltsFace],
                                       slipRateMagnitude[ltsFace],
                                       startTime);

        static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace, startTime);

        common::savePeakSlipRateOutput(slipRateMagnitude[ltsFace], peakSlipRate[ltsFace]);

        if (this->drParameters->isFrictionEnergyRequired &&
            this->drParameters->isCheckAbortCriteraEnabled) {
          common::updateTimeSinceSlipRateBelowThreshold(
              slipRateMagnitude[ltsFace],
              ruptureTimePending[ltsFace],
              energyData[ltsFace],
              this->deltaT[timeIndex],
              this->drParameters->terminatorSlipRateThreshold);
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
                                                   impAndEta[ltsFace],
                                                   impedanceMatrices[ltsFace],
                                                   imposedStatePlus[ltsFace],
                                                   imposedStateMinus[ltsFace],
                                                   qInterpolatedPlus[ltsFace],
                                                   qInterpolatedMinus[ltsFace],
                                                   timeWeights);
      LIKWID_MARKER_STOP("computeDynamicRupturePostcomputeImposedState");
      SCOREP_USER_REGION_END(myRegionHandle)

      if (this->drParameters->isFrictionEnergyRequired) {
        common::computeFrictionEnergy(energyData[ltsFace],
                                      qInterpolatedPlus[ltsFace],
                                      qInterpolatedMinus[ltsFace],
                                      impAndEta[ltsFace],
                                      timeWeights,
                                      spaceWeights,
                                      godunovData[ltsFace],
                                      slipRateMagnitude[ltsFace],
                                      this->drParameters->energiesFromAcrossFaultVelocities);
      }
    }
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
