// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_

#include <Common/Constants.h>
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
  public:
  explicit BaseFrictionLaw(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolver(drParameters) {};

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializer::Layer& layerData,
                const seissol::initializer::DynamicRupture* const dynRup,
                real fullUpdateTime,
                const double timeWeights[ConvergenceOrder],
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (layerData.getNumberOfCells() == 0) {
      return;
    }
    auto& self = *this;
    runtime.enqueueHost([fullUpdateTime, dynRup, &self, &layerData] {
      static_cast<BaseFrictionLaw&>(self).copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
      static_cast<Derived&>(self).copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    });
    runtime.enqueueOmpFor(
        layerData.getNumberOfCells(),
        [fullUpdateTime, &self, &layerData, timeWeights](std::size_t ltsFace) {
          alignas(Alignment) FaultStresses<Executor::Host> faultStresses{};
          common::precomputeStressFromQInterpolated(faultStresses,
                                                    self.impAndEta[ltsFace],
                                                    self.impedanceMatrices[ltsFace],
                                                    self.qInterpolatedPlus[ltsFace],
                                                    self.qInterpolatedMinus[ltsFace]);

          // define some temporary variables
          std::array<real, misc::NumPaddedPoints> stateVariableBuffer{0};
          std::array<real, misc::NumPaddedPoints> strengthBuffer{0};

          static_cast<Derived&>(self).preHook(stateVariableBuffer, ltsFace);
          TractionResults<Executor::Host> tractionResults = {};

          // loop over sub time steps (i.e. quadrature points in time)
          real updateTime = self.mFullUpdateTime;
          for (std::size_t timeIndex = 0; timeIndex < ConvergenceOrder; timeIndex++) {
            updateTime += self.deltaT[timeIndex];
            common::adjustInitialStress(self.initialStressInFaultCS[ltsFace],
                                        self.nucleationStressInFaultCS[ltsFace],
                                        self.initialPressure[ltsFace],
                                        self.nucleationPressure[ltsFace],
                                        updateTime,
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

          common::savePeakSlipRateOutput(self.slipRateMagnitude[ltsFace],
                                         self.peakSlipRate[ltsFace]);
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
                                          self.spaceWeights,
                                          self.godunovData[ltsFace],
                                          self.slipRateMagnitude[ltsFace],
                                          self.drParameters->energiesFromAcrossFaultVelocities);
          }
        });
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_BASEFRICTIONLAW_H_
