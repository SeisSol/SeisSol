// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_

#include "DynamicRupture/Output/ReceiverBasedOutput.h"

namespace seissol::dr::output {
struct RateAndState : public ReceiverOutputImpl<RateAndState> {
  template<typename Cfg>
  Real<Cfg> computeLocalStrength(LocalInfo<Cfg>& local) {
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
           std::min(effectiveNormalStress, static_cast<Real<Cfg>>(0.0));
  }

  template<typename Cfg>
  Real<Cfg> computeStateVariable(LocalInfo<Cfg>& local) {
    return getCellData<LTSRateAndState::StateVariable>(Cfg(), local)[local.gpIndex];
  }

  template<typename Cfg>
  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const {
    auto baseVector = ReceiverOutputImpl::getOutputVariables();
    baseVector.push_back(drStorage->info<LTSRateAndState::StateVariable>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
