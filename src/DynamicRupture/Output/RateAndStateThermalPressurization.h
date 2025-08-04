// SPDX-FileCopyrightText: 2021 SeisSol Group
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
struct RateAndStateThermalPressurization
    : public ReceiverOutputImpl<RateAndStateThermalPressurization> {
  public:
  template <typename Cfg>
  Real<Cfg> computeLocalStrength(LocalInfo<Cfg>& local) {
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
           std::min(effectiveNormalStress, static_cast<Real<Cfg>>(0.0));
  }

  template <typename Cfg>
  Real<Cfg> computeStateVariable(LocalInfo<Cfg>& local) {
    return getCellData<LTSRateAndState::StateVariable>(Cfg(), local)[local.gpIndex];
  }

  template <typename Cfg>
  Real<Cfg> computeFluidPressure(LocalInfo<Cfg>& local) {
    const auto* const pressure =
        ReceiverOutputImpl::getCellData<LTSThermalPressurization::Pressure>(Cfg(), local);
    return pressure[local.gpIndex];
  }

  template <typename Cfg>
  void outputSpecifics(const std::shared_ptr<ReceiverOutputData>& outputData,
                       const LocalInfo<Cfg>& local,
                       size_t cacheLevel,
                       size_t receiverIdx) {
    auto& tpVariables = std::get<VariableID::ThermalPressurizationVariables>(outputData->vars);
    if (tpVariables.isActive) {
      const auto* const temperature =
          getCellData<LTSThermalPressurization::Temperature>(Cfg(), local);
      tpVariables(TPID::Temperature, cacheLevel, receiverIdx) = temperature[local.gpIndex];

      const auto* const pressure = getCellData<LTSThermalPressurization::Pressure>(Cfg(), local);
      tpVariables(TPID::Pressure, cacheLevel, receiverIdx) = pressure[local.gpIndex];
    }
  }

  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = ReceiverOutputImpl::getOutputVariables();
    // manually add the R+S state var here for now
    baseVector.push_back(drStorage->info<LTSRateAndState::StateVariable>().index);
    baseVector.push_back(this->drStorage->info<LTSThermalPressurization::Temperature>().index);
    baseVector.push_back(this->drStorage->info<LTSThermalPressurization::Pressure>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_RATEANDSTATETHERMALPRESSURIZATION_H_
