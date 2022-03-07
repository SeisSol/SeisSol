#ifndef SEISSOL_DR_OUTPUT_RS_TP_HPP
#define SEISSOL_DR_OUTPUT_RS_TP_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class RateAndStateThermalPressurisation : public RateAndState {
  public:
  using RateAndState::postCompute;
  using RateAndState::RateAndState;

  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& eInteroperability) override {
    RateAndState::tiePointers(layerData, dynRup, eInteroperability);

    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);

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
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_RS_TP_HPP
