#ifndef SEISSOL_DR_OUTPUT_RS_TP_HPP
#define SEISSOL_DR_OUTPUT_RS_TP_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class RateAndStateThermalPressurization : public RateAndState {
  public:
  using RateAndState::RateAndState;

  protected:
  real computeFluidPressure(LocalInfo& local) override {
    using DrLtsDescrType = seissol::initializer::LTSRateAndStateThermalPressurization;
    auto const* const pressure = local.layer->var(static_cast<DrLtsDescrType*>(drDescr)->pressure);
    return pressure[local.ltsId][local.nearestGpIndex];
  }
  void outputSpecifics(std::shared_ptr<ReceiverOutputData>& outputData,
                       const LocalInfo& local,
                       size_t cacheLevel,
                       size_t receiverIdx) override {
    auto& tpVariables = std::get<VariableID::ThermalPressurizationVariables>(outputData->vars);
    if (tpVariables.isActive) {
      using DrLtsDescrType = seissol::initializer::LTSRateAndStateThermalPressurization;
      auto const* const temperature =
          local.layer->var(static_cast<DrLtsDescrType*>(drDescr)->temperature);
      tpVariables(TPID::Temperature, cacheLevel, receiverIdx) =
          temperature[local.ltsId][local.nearestGpIndex];

      auto const* const pressure =
          local.layer->var(static_cast<DrLtsDescrType*>(drDescr)->pressure);
      tpVariables(TPID::Pressure, cacheLevel, receiverIdx) =
          pressure[local.ltsId][local.nearestGpIndex];
    }
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_RS_TP_HPP
