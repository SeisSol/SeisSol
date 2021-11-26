#ifndef SEISSOL_DR_OUTPUT_RS_HPP
#define SEISSOL_DR_OUTPUT_RS_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class RateAndState : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& e_interoperability) override {
    Base::tiePointers(layerData, dynRup, e_interoperability);

    auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto numGaussPoints2d = init::QInterpolated::Stop[0];
    real(*slipRateStrike)[numGaussPoints2d] = layerData.var(concreteLts->slipRateStrike);
    real(*slipRateDip)[numGaussPoints2d] = layerData.var(concreteLts->slipRateDip);
    real(*mu)[numGaussPoints2d] = layerData.var(concreteLts->mu);
    real(*stateVar)[numGaussPoints2d] = layerData.var(concreteLts->stateVariable);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRateStrike, slipRateDip, mu);
      e_interoperability.copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVar);
    }
  }

  void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_RS_HPP
