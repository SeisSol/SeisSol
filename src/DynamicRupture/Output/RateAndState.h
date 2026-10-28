// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_

#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::output {
/**
  The output of the rate-and-state friction laws, for `Derived`, which may record more (CRTP).
 */
template <typename Derived>
class RateAndStateImpl : public ReceiverOutputImpl<Derived> {
  public:
  template <typename Cfg>
  using LocalInfo = typename ReceiverOutputImpl<Derived>::template LocalInfo<Cfg>;

  template <typename Cfg>
  Real<Cfg> computeLocalStrength(LocalInfo<Cfg>& local) {
    using real = Real<Cfg>;
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
           std::min(effectiveNormalStress, static_cast<real>(0.0));
  }

  template <typename Cfg>
  Real<Cfg> computeLocalStrengthSlope(LocalInfo<Cfg>& local) {
    using real = Real<Cfg>;
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return effectiveNormalStress < 0 ? local.frictionCoefficient : static_cast<real>(0.0);
  }

  template <typename Cfg>
  Real<Cfg> computeStateVariable(LocalInfo<Cfg>& local) {
    return this->template getCellData<LTSRateAndState::StateVariable>(local)[local.gpIndex];
  }

  template <typename Cfg>
  void handleNonConvergence(LocalInfo<Cfg>& local) {
    const auto* inner = this->template getCellData<LTSRateAndState::ConvergenceInner>(local);
    const auto* outer = this->template getCellData<LTSRateAndState::ConvergenceOuter>(local);
    std::vector<std::size_t> failuresInner;
    std::vector<std::size_t> failuresOuter;
    for (std::size_t i = 0; i < misc::NumBoundaryGaussPoints<Cfg>; ++i) {
      const auto index = i * Cfg::NumSimulations + local.fusedIndex;
      if (!inner[index]) {
        failuresInner.push_back(i);
      }
      if (!outer[index]) {
        failuresOuter.push_back(i);
      }
    }

    if (!(failuresInner.empty() && failuresOuter.empty())) {
      const auto& point = local.state->receivers[local.index].global;
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

class RateAndState : public RateAndStateImpl<RateAndState> {};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
