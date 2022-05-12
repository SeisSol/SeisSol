#ifndef SEISSOL_DR_OUTPUT_RS_TP_HPP
#define SEISSOL_DR_OUTPUT_RS_TP_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class RateAndStateThermalPressurization : public RateAndState {
  public:
  using RateAndState::RateAndState;
  using RateAndState::tiePointers;

  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& eInteroperability) override {
    RateAndState::tiePointers(layerData, dynRup, eInteroperability);

    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurization*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real(*fluidPressure)[misc::numPaddedPoints] = layerData.var(concreteLts->pressure);
    real(*fluidTemperature)[misc::numPaddedPoints] = layerData.var(concreteLts->temperature);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      eInteroperability.copyFrictionOutputToFortranThermalPressurization(
          ltsFace, meshFace, fluidPressure, fluidTemperature);
    }
  }

  protected:
  real computePf() override { return 0.0; }
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
