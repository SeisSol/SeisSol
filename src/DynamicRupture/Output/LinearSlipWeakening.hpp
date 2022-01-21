#ifndef SEISSOL_DR_OUTPUT_LSW_HPP
#define SEISSOL_DR_OUTPUT_LSW_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class LinearSlipWeakening : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& eInteroperability) override {
    Base::tiePointers(layerData, dynRup, eInteroperability);

    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*slipRate1)[size] = layerData.var(concreteLts->slipRate1);
    real(*slipRate2)[size] = layerData.var(concreteLts->slipRate2);
    real(*mu)[size] = layerData.var(concreteLts->mu);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      eInteroperability.copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRate1, slipRate2, mu);
    }
  }

  void postCompute(seissol::initializers::DynamicRupture& dynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_LSW_HPP
