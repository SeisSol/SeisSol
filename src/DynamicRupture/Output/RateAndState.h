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

  std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = ReceiverOutput::getOutputVariables();
    baseVector.push_back(drStorage->info<LTSRateAndState::StateVariable>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
