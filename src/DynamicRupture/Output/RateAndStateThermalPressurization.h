// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATETHERMALPRESSURIZATION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATETHERMALPRESSURIZATION_H_

#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include <DynamicRupture/Output/RateAndState.h>
#include <Memory/Descriptor/DynamicRupture.h>

namespace seissol::dr::output {
class RateAndStateThermalPressurization : public RateAndState {
  public:
  using RateAndState::RateAndState;

  protected:
  real computeFluidPressure(LocalInfo& local) override {
    using DrLtsDescrType = seissol::initializer::ThermalPressurization;
    const auto* const pressure =
        getCellData(local, dynamic_cast<DrLtsDescrType*>(drDescr)->pressure);
    return pressure[local.nearestGpIndex];
  }
  void outputSpecifics(std::shared_ptr<ReceiverOutputData>& outputData,
                       const LocalInfo& local,
                       size_t cacheLevel,
                       size_t receiverIdx) override {
    auto& tpVariables = std::get<VariableID::ThermalPressurizationVariables>(outputData->vars);
    if (tpVariables.isActive) {
      using DrLtsDescrType = seissol::initializer::ThermalPressurization;
      const auto* const temperature =
          getCellData(local, dynamic_cast<DrLtsDescrType*>(drDescr)->temperature);
      tpVariables(TPID::Temperature, cacheLevel, receiverIdx) = temperature[local.nearestGpIndex];

      const auto* const pressure =
          getCellData(local, dynamic_cast<DrLtsDescrType*>(drDescr)->pressure);
      tpVariables(TPID::Pressure, cacheLevel, receiverIdx) = pressure[local.nearestGpIndex];
    }
  }

  std::vector<std::size_t> getOutputVariables() const override {
    using DrLtsDescrType = seissol::initializer::ThermalPressurization;
    auto baseVector = RateAndState::getOutputVariables();
    baseVector.push_back(dynamic_cast<const DrLtsDescrType*>(drDescr)->temperature.index);
    baseVector.push_back(dynamic_cast<const DrLtsDescrType*>(drDescr)->pressure.index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATETHERMALPRESSURIZATION_H_
