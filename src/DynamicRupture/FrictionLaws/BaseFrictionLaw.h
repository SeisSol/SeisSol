#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "FrictionSolver.h"
#include "FrictionSolverCommon.h"

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
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      alignas(ALIGNMENT) FaultStresses faultStresses{};
      LIKWID_MARKER_START("PrecomputeStress");
      common::precomputeStressFromQInterpolated(faultStresses,
                                                impAndEta[ltsFace],
                                                qInterpolatedPlus[ltsFace],
                                                qInterpolatedMinus[ltsFace]);
      LIKWID_MARKER_STOP("PrecomputeStress");

      LIKWID_MARKER_START("PreHook");
      // define some temporary variables
      std::array<real, misc::numPaddedPoints> stateVariableBuffer{0};
      std::array<real, misc::numPaddedPoints> strengthBuffer{0};

      static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);
      LIKWID_MARKER_STOP("PreHook");

      LIKWID_MARKER_START("UpdateFrictionAndSlip");
      TractionResults tractionResults = {};

      // loop over sub time steps (i.e. quadrature points in time)
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
        this->adjustInitialStress(ltsFace, timeIndex);
        static_cast<Derived*>(this)->updateFrictionAndSlip(faultStresses,
                                                           tractionResults,
                                                           stateVariableBuffer,
                                                           strengthBuffer,
                                                           ltsFace,
                                                           timeIndex);
      }
      LIKWID_MARKER_STOP("UpdateFrictionAndSlip");

      LIKWID_MARKER_START("PostHook");
      static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);

      common::saveRuptureFrontOutput(ruptureTimePending[ltsFace],
                                     ruptureTime[ltsFace],
                                     slipRateMagnitude[ltsFace],
                                     mFullUpdateTime);

      static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace);

      common::savePeakSlipRateOutput(slipRateMagnitude[ltsFace], peakSlipRate[ltsFace]);
      LIKWID_MARKER_STOP("PostHook");

      LIKWID_MARKER_START("PostcomputeImposedState");
      common::postcomputeImposedStateFromNewStress(faultStresses,
                                                   tractionResults,
                                                   impAndEta[ltsFace],
                                                   imposedStatePlus[ltsFace],
                                                   imposedStateMinus[ltsFace],
                                                   qInterpolatedPlus[ltsFace],
                                                   qInterpolatedMinus[ltsFace],
                                                   timeWeights);
      LIKWID_MARKER_STOP("PostcomputeImposedState");
    }
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_BASEFRICTIONLAW_H
