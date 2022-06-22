#ifndef SEISSOL_DR_OUTPUT_RS_TP_HPP
#define SEISSOL_DR_OUTPUT_RS_TP_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class RateAndStateThermalPressurization : public RateAndState {
  public:
  using RateAndState::RateAndState;

  protected:
  real computeFluidPressure() override {
    using DrLtsDescrT = seissol::initializers::LTS_RateAndStateThermalPressurization;
    auto* pressure = local.layer->var(static_cast<DrLtsDescrT*>(drDescr)->pressure);
    return pressure[local.ltsId][local.nearestGpIndex];
  }
  void outputSpecifics(ReceiverBasedOutputData& outputData,
                       size_t cacheLevel,
                       size_t receiverIdx) override {
    auto& tpVariables = std::get<VariableID::TpVariables>(outputData.vars);
    if (tpVariables.isActive) {
      using DrLtsDescrT = seissol::initializers::LTS_RateAndStateThermalPressurization;
      auto* temperature = local.layer->var(static_cast<DrLtsDescrT*>(drDescr)->temperature);
      tpVariables(TPID::Temperature, cacheLevel, receiverIdx) =
          temperature[local.ltsId][local.nearestGpIndex];

      auto* pressure = local.layer->var(static_cast<DrLtsDescrT*>(drDescr)->pressure);
      tpVariables(TPID::Pressure, cacheLevel, receiverIdx) =
          pressure[local.ltsId][local.nearestGpIndex];
    }
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_RS_TP_HPP
