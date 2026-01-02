// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATETHERMALPRESSURIZATION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATETHERMALPRESSURIZATION_H_

#include "DynamicRupture/Output/RateAndState.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::output {
class RateAndStateThermalPressurization : public RateAndState {
  public:
  using RateAndState::RateAndState;

  protected:
  real computeFluidPressure(LocalInfo& local) override {
    const auto* const pressure = getCellData<LTSThermalPressurization::Pressure>(local);
    return pressure[local.gpIndex];
  }
  void outputSpecifics(const std::shared_ptr<ReceiverOutputData>& outputData,
                       const LocalInfo& local,
                       size_t cacheLevel,
                       size_t receiverIdx) override {
    auto& tpVariables = std::get<VariableID::ThermalPressurizationVariables>(outputData->vars);
    if (tpVariables.isActive) {
      const auto* const temperature = getCellData<LTSThermalPressurization::Temperature>(local);
      tpVariables(TPID::Temperature, cacheLevel, receiverIdx) = temperature[local.gpIndex];

      const auto* const pressure = getCellData<LTSThermalPressurization::Pressure>(local);
      tpVariables(TPID::Pressure, cacheLevel, receiverIdx) = pressure[local.gpIndex];
    }
  }

  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = RateAndState::getOutputVariables();
    baseVector.push_back(drStorage->info<LTSThermalPressurization::Temperature>().index);
    baseVector.push_back(drStorage->info<LTSThermalPressurization::Pressure>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATETHERMALPRESSURIZATION_H_
