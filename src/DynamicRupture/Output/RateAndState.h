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
  std::array<real, seissol::multisim::NumSimulations> computeLocalStrength(LocalInfo& local) override {
    std::array<real, seissol::multisim::NumSimulations> strength{};
    for (unsigned int i = 0; i < seissol::multisim::NumSimulations; ++i) {
      const auto effectiveNormalStress =
        local.transientNormalTraction[i] + local.iniNormalTraction[i] - local.fluidPressure[i];
      strength[i] = -1.0 * local.frictionCoefficient[i] *
                     std::min(effectiveNormalStress, static_cast<real>(0.0));
      }
      return strength;
  }

  real computeStateVariable(LocalInfo& local, unsigned int index) override {
    const auto* descr = reinterpret_cast<seissol::initializer::LTSRateAndState*>(drDescr);
    assert((descr != nullptr) && "dr descr. must be a subtype of LTS_RateAndState");
    return getCellData(local, descr->stateVariable)[index];
  }

  std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = ReceiverOutput::getOutputVariables();
    baseVector.push_back(
        static_cast<seissol::initializer::LTSRateAndState*>(drDescr)->stateVariable.index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATE_H_
