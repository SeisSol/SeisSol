#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "FrictionSolver.h"
#include "FrictionSolverCommon.h"
#include "Monitoring/instrumentation.hpp"
#include "Common/constants.hpp"

namespace seissol::dr::friction_law {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Config, typename Derived>
class BaseFrictionLaw : public FrictionSolver<Config> {
  public:
  using RealT = typename Config::RealT;
  explicit BaseFrictionLaw(dr::DRParameters* drParameters) : FrictionSolver<Config>(drParameters){};

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture<Config> const* const dynRup,
                RealT fullUpdateTime,
                const double timeWeights[Config::ConvergenceOrder]) override {
    SCOREP_USER_REGION_DEFINE(myRegionHandle)
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      alignas(Alignment) FaultStresses<Config> faultStresses{};
      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePrecomputeStress", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePrecomputeStress");
      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace]);
      LIKWID_MARKER_STOP("computeDynamicRupturePrecomputeStress");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePreHook", SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePreHook");
      // define some temporary variables
      std::array<RealT, misc::numPaddedPoints<Config>> stateVariableBuffer{0};
      std::array<RealT, misc::numPaddedPoints<Config>> strengthBuffer{0};

      static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
      LIKWID_MARKER_STOP("computeDynamicRupturePreHook");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRuptureUpdateFrictionAndSlip",
                               SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRuptureUpdateFrictionAndSlip");
      TractionResults<Config> tractionResults = {};

      // loop over sub time steps (i.e. quadrature points in time)
      for (unsigned timeIndex = 0; timeIndex < Config::ConvergenceOrder; timeIndex++) {
        common::adjustInitialStress(initialStressInFaultCS[ltsFace],
                                    nucleationStressInFaultCS[ltsFace],
                                    this->mFullUpdateTime,
                                    this->drParameters->t0,
                                    this->deltaT[timeIndex]);

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

      common::saveRuptureFrontOutput(ruptureTimePending[ltsFace],
                                     ruptureTime[ltsFace],
                                     slipRateMagnitude[ltsFace],
                                     mFullUpdateTime);

      static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace);

      common::savePeakSlipRateOutput(slipRateMagnitude[ltsFace], peakSlipRate[ltsFace]);
      LIKWID_MARKER_STOP("computeDynamicRupturePostHook");
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRupturePostcomputeImposedState",
                               SCOREP_USER_REGION_TYPE_COMMON)
      LIKWID_MARKER_START("computeDynamicRupturePostcomputeImposedState");
      common::postcomputeImposedStateFromNewStress(faultStresses,
                                                   tractionResults,
                                                   impAndEta[ltsFace],
                                                   imposedStatePlus[ltsFace],
                                                   imposedStateMinus[ltsFace],
                                                   qInterpolatedPlus[ltsFace],
                                                   qInterpolatedMinus[ltsFace],
                                                   timeWeights);
      LIKWID_MARKER_STOP("computeDynamicRupturePostcomputeImposedState");
      SCOREP_USER_REGION_END(myRegionHandle)

      common::computeFrictionEnergy(energyData[ltsFace],
                                    qInterpolatedPlus[ltsFace],
                                    qInterpolatedMinus[ltsFace],
                                    impAndEta[ltsFace],
                                    timeWeights,
                                    spaceWeights,
                                    godunovData[ltsFace]);
    }
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_BASEFRICTIONLAW_H
