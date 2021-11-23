#ifndef SEISSOL_DR_OUTOUT_LSW_HPP
#define SEISSOL_DR_OUTOUT_LSW_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class LinearSlipWeakening : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& e_interoperability) override {
    Base::tiePointers(layerData, dynRup, e_interoperability);

    auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*slipRateStrike)[size] = layerData.var(concreteLts->slipRateStrike);
    real(*slipRateDip)[size] = layerData.var(concreteLts->slipRateDip);
    real(*mu)[size] = layerData.var(concreteLts->mu);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRateStrike, slipRateDip, mu);
    }
  }

  void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTOUT_LSW_HPP
