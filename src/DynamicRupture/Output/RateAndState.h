// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::output {
class RateAndState : public ReceiverOutput {
  protected:
  real computeLocalStrength(LocalInfo& local) override {
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
           std::min(effectiveNormalStress, static_cast<real>(0.0));
  }

  real computeStateVariable(LocalInfo& local) override {
    return getCellData<LTSRateAndState::StateVariable>(local)[local.gpIndex];
  }

  void handleNonConvergence(LocalInfo& local) override {
    const auto* inner = getCellData<LTSRateAndState::ConvergenceInner>(local);
    const auto* outer = getCellData<LTSRateAndState::ConvergenceOuter>(local);
    std::vector<std::size_t> failuresInner;
    std::vector<std::size_t> failuresOuter;
    for (std::size_t i = 0; i < misc::NumBoundaryGaussPoints; ++i) {
      if (!inner[i]) {
        failuresInner.push_back(i);
      }
      if (!outer[i]) {
        failuresOuter.push_back(i);
      }
    }

    if (!(failuresInner.empty() && failuresOuter.empty())) {
      const auto* pointData = local.state->receiverPoints[local.index].global.coords;
      const std::array<double, 3> point{pointData[0], pointData[1], pointData[2]};
      logError() << "A rate and state cell failed to converge; at the cell around" << point
                 << ". PointIDs of failure (inner, outer loop failures):" << failuresInner
                 << failuresOuter;
    }
  }

  public:
  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = ReceiverOutput::getOutputVariables();
    baseVector.push_back(drStorage->info<LTSRateAndState::StateVariable>().index);
    baseVector.push_back(drStorage->info<LTSRateAndState::ConvergenceInner>().index);
    baseVector.push_back(drStorage->info<LTSRateAndState::ConvergenceOuter>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
