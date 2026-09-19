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
/**
  The parts which rate-and-state shares with its thermal-pressurization variant. Kept as a mixin
  rather than a base class of the latter, so that the customisation points still resolve to the
  most derived implementation.
 */
template <typename Derived>
class RateAndStateBase : public ReceiverOutputImpl<Derived> {
  public:
  using LocalInfo = ReceiverOutput::LocalInfo;

  real computeLocalStrength(LocalInfo& local) {
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
           std::min(effectiveNormalStress, static_cast<real>(0.0));
  }

  real computeStateVariable(LocalInfo& local) {
    return this->template getCellData<LTSRateAndState::StateVariable>(local)[local.gpIndex];
  }

  void handleNonConvergence(LocalInfo& local) {
    const auto* inner = this->template getCellData<LTSRateAndState::ConvergenceInner>(local);
    const auto* outer = this->template getCellData<LTSRateAndState::ConvergenceOuter>(local);
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
      const auto* pointData = local.state->receivers[local.index].global.coords;
      const std::array<double, 3> point{pointData[0], pointData[1], pointData[2]};
      auto& printWarning = *local.printWarning;
      if (!printWarning) {
        logWarning(true)
            << "A rate and state cell failed to converge at the given settings; at the cell around"
            << point << "at simulation time" << local.time
            << "s. PointIDs of failure (inner, outer loop failures):" << failuresInner
            << failuresOuter;
        printWarning = true;
      }
    }
  }

  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = ReceiverOutput::getOutputVariables();
    baseVector.push_back(this->drStorage_->template info<LTSRateAndState::StateVariable>().index);
    baseVector.push_back(
        this->drStorage_->template info<LTSRateAndState::ConvergenceInner>().index);
    baseVector.push_back(
        this->drStorage_->template info<LTSRateAndState::ConvergenceOuter>().index);
    return baseVector;
  }
};

class RateAndState : public RateAndStateBase<RateAndState> {};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
