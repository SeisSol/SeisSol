#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "FrictionSolver.h"
#include "FrictionSolverCommon.h"
#include "Monitoring/instrumentation.fpp"

namespace seissol::dr::friction_law {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  public:
  explicit BaseFrictionLaw(dr::DRParameters* drParameters) : FrictionSolver(drParameters){};

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture const* const dynRup,
                real fullUpdateTime,
                const double timeWeights[CONVERGENCE_ORDER]) override {
    SCOREP_USER_REGION_DEFINE(myRegionHandle)
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      alignas(ALIGNMENT) FaultStresses faultStresses{};
      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePrecomputeStress", SCOREP_USER_REGION_TYPE_COMMON)
      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace]);
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePreHook", SCOREP_USER_REGION_TYPE_COMMON)
      // define some temporary variables
      std::array<real, misc::numPaddedPoints> stateVariableBuffer{0};
      std::array<real, misc::numPaddedPoints> strengthBuffer{0};

      static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRuptureUpdateFrictionAndSlip",
                               SCOREP_USER_REGION_TYPE_COMMON)
      TractionResults tractionResults = {};

      // loop over sub time steps (i.e. quadrature points in time)
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
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
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(
          myRegionHandle, "computeDynamicRupturePostHook", SCOREP_USER_REGION_TYPE_COMMON)
      static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);

      common::saveRuptureFrontOutput(ruptureTimePending[ltsFace],
                                     ruptureTime[ltsFace],
                                     slipRateMagnitude[ltsFace],
                                     mFullUpdateTime);

      static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace);

      common::savePeakSlipRateOutput(slipRateMagnitude[ltsFace], peakSlipRate[ltsFace]);
      SCOREP_USER_REGION_END(myRegionHandle)

      SCOREP_USER_REGION_BEGIN(myRegionHandle,
                               "computeDynamicRupturePostcomputeImposedState",
                               SCOREP_USER_REGION_TYPE_COMMON)
      common::postcomputeImposedStateFromNewStress(faultStresses,
                                                   tractionResults,
                                                   impAndEta[ltsFace],
                                                   imposedStatePlus[ltsFace],
                                                   imposedStateMinus[ltsFace],
                                                   qInterpolatedPlus[ltsFace],
                                                   qInterpolatedMinus[ltsFace],
                                                   timeWeights);
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
