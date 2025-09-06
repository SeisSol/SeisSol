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
template <typename Cfg, typename Derived>
class BaseFrictionLaw : public FrictionSolverImpl<Cfg> {
  private:
  size_t currLayerSize{};

  public:
  using real = Real<Cfg>;
  explicit BaseFrictionLaw(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolverImpl<Cfg>(drParameters) {}

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
  void evaluate(double fullUpdateTime,
                const FrictionSolver::FrictionTime& frictionTime,
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
      alignas(Alignment) FaultStresses<Cfg, Executor::Host> faultStresses{};
      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePrecomputeStress", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePrecomputeStress");
      const auto etaPDamp =
          this->drParameters->etaStop > this->mFullUpdateTime ? this->drParameters->etaHack : 1.0;
      common::precomputeStressFromQInterpolated<Cfg>(faultStresses,
                                                     this->impAndEta[ltsFace],
                                                     this->impedanceMatrices[ltsFace],
                                                     this->qInterpolatedPlus[ltsFace],
                                                     this->qInterpolatedMinus[ltsFace],
                                                     etaPDamp);
      LIKWID_MARKER_STOP("computeDynamicRupturePrecomputeStress");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePreHook", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePreHook");
      // define some temporary variables
      std::array<real, misc::NumPaddedPoints<Cfg>> stateVariableBuffer{0};
      std::array<real, misc::NumPaddedPoints<Cfg>> strengthBuffer{0};

      static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
      LIKWID_MARKER_STOP("computeDynamicRupturePreHook");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRuptureUpdateFrictionAndSlip",
                               SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRuptureUpdateFrictionAndSlip");
      TractionResults<Cfg, Executor::Host> tractionResults = {};

      // loop over sub time steps (i.e. quadrature points in time)
      real updateTime = this->mFullUpdateTime;
      for (std::size_t timeIndex = 0; timeIndex < Cfg::ConvergenceOrder; timeIndex++) {
        updateTime += this->deltaT[timeIndex];
        for (unsigned i = 0; i < this->drParameters->nucleationCount; ++i) {
          common::adjustInitialStress<Cfg>(
              this->initialStressInFaultCS[ltsFace],
              this->nucleationStressInFaultCS[ltsFace * this->drParameters->nucleationCount + i],
              this->initialPressure[ltsFace],
              this->nucleationPressure[ltsFace * this->drParameters->nucleationCount + i],
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
      }
      LIKWID_MARKER_STOP("computeDynamicRuptureUpdateFrictionAndSlip");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePostHook", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePostHook");
      static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);

      common::saveRuptureFrontOutput<Cfg>(this->ruptureTimePending[ltsFace],
                                          this->ruptureTime[ltsFace],
                                          this->slipRateMagnitude[ltsFace],
                                          this->mFullUpdateTime);

      static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace);

      common::savePeakSlipRateOutput<Cfg>(this->slipRateMagnitude[ltsFace],
                                          this->peakSlipRate[ltsFace]);
      LIKWID_MARKER_STOP("computeDynamicRupturePostHook");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRupturePostcomputeImposedState",
                               SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePostcomputeImposedState");
      common::postcomputeImposedStateFromNewStress<Cfg>(faultStresses,
                                                        tractionResults,
                                                        this->impAndEta[ltsFace],
                                                        this->impedanceMatrices[ltsFace],
                                                        this->imposedStatePlus[ltsFace],
                                                        this->imposedStateMinus[ltsFace],
                                                        this->qInterpolatedPlus[ltsFace],
                                                        this->qInterpolatedMinus[ltsFace],
                                                        timeWeights);
      LIKWID_MARKER_STOP("computeDynamicRupturePostcomputeImposedState");
      SCOREP_USER_REGION_END(myRegionHandle)

      if (this->drParameters->isFrictionEnergyRequired) {

        if (this->drParameters->isCheckAbortCriteraEnabled) {
          common::updateTimeSinceSlipRateBelowThreshold<Cfg>(
              this->slipRateMagnitude[ltsFace],
              this->ruptureTimePending[ltsFace],
              this->energyData[ltsFace],
              this->sumDt,
              this->drParameters->terminatorSlipRateThreshold);
        }
        common::computeFrictionEnergy<Cfg>(this->energyData[ltsFace],
                                           this->qInterpolatedPlus[ltsFace],
                                           this->qInterpolatedMinus[ltsFace],
                                           this->impAndEta[ltsFace],
                                           timeWeights,
                                           this->spaceWeights,
                                           this->godunovData[ltsFace],
                                           this->slipRateMagnitude[ltsFace],
                                           this->drParameters->energiesFromAcrossFaultVelocities);
      }
    }
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
