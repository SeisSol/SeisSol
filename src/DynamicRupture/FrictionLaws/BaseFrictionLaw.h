#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Misc.h"
#include "FrictionSolver.h"
#include "FrictionSolverCommon.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Monitoring/instrumentation.hpp"

namespace seissol::dr::friction_law {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  public:
  explicit BaseFrictionLaw(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolver(drParameters) {};

  void dependency(seissol::parallel::runtime::StreamRuntime& runtime) override {
    // TODO
  }

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializer::Layer& layerData,
                const seissol::initializer::DynamicRupture* const dynRup,
                real fullUpdateTime,
                const double timeWeights[CONVERGENCE_ORDER],
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    auto& self = *this;
    runtime.enqueueHost([=, &self, &layerData] {
      BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
      static_cast<Derived&>(self).copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

      // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
        alignas(ALIGNMENT) FaultStresses faultStresses{};
        common::precomputeStressFromQInterpolated(faultStresses,
                                                  self.impAndEta[ltsFace],
                                                  self.impedanceMatrices[ltsFace],
                                                  self.qInterpolatedPlus[ltsFace],
                                                  self.qInterpolatedMinus[ltsFace]);

        // define some temporary variables
        std::array<real, misc::numPaddedPoints> stateVariableBuffer{0};
        std::array<real, misc::numPaddedPoints> strengthBuffer{0};

        static_cast<Derived&>(self).preHook(stateVariableBuffer, ltsFace);
        TractionResults tractionResults = {};

        // loop over sub time steps (i.e. quadrature points in time)
        for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
          common::adjustInitialStress(self.initialStressInFaultCS[ltsFace],
                                      self.nucleationStressInFaultCS[ltsFace],
                                      self.initialPressure[ltsFace],
                                      self.nucleationPressure[ltsFace],
                                      self.mFullUpdateTime,
                                      self.drParameters->t0,
                                      self.deltaT[timeIndex]);

          static_cast<Derived&>(self).updateFrictionAndSlip(faultStresses,
                                                            tractionResults,
                                                            stateVariableBuffer,
                                                            strengthBuffer,
                                                            ltsFace,
                                                            timeIndex);
        }
        static_cast<Derived&>(self).postHook(stateVariableBuffer, ltsFace);

        common::saveRuptureFrontOutput(self.ruptureTimePending[ltsFace],
                                       self.ruptureTime[ltsFace],
                                       self.slipRateMagnitude[ltsFace],
                                       self.mFullUpdateTime);

        static_cast<Derived&>(self).saveDynamicStressOutput(ltsFace);

        common::savePeakSlipRateOutput(self.slipRateMagnitude[ltsFace], self.peakSlipRate[ltsFace]);
        common::postcomputeImposedStateFromNewStress(faultStresses,
                                                     tractionResults,
                                                     self.impAndEta[ltsFace],
                                                     self.impedanceMatrices[ltsFace],
                                                     self.imposedStatePlus[ltsFace],
                                                     self.imposedStateMinus[ltsFace],
                                                     self.qInterpolatedPlus[ltsFace],
                                                     self.qInterpolatedMinus[ltsFace],
                                                     timeWeights);

        if (self.drParameters->isFrictionEnergyRequired) {

          if (self.drParameters->isCheckAbortCriteraEnabled) {
            common::updateTimeSinceSlipRateBelowThreshold(
                self.slipRateMagnitude[ltsFace],
                self.ruptureTimePending[ltsFace],
                self.energyData[ltsFace],
                self.sumDt,
                self.drParameters->terminatorSlipRateThreshold);
          }
          common::computeFrictionEnergy(self.energyData[ltsFace],
                                        self.qInterpolatedPlus[ltsFace],
                                        self.qInterpolatedMinus[ltsFace],
                                        self.impAndEta[ltsFace],
                                        timeWeights,
                                        spaceWeights,
                                        self.godunovData[ltsFace]);
        }
      }
    });
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_BASEFRICTIONLAW_H
