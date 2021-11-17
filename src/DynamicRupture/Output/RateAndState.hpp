#ifndef SEISSOL_DR_OUTOUT_RS_AGING_LAW_HPP
#define SEISSOL_DR_OUTOUT_RS_AGING_LAW_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class RateAndState : public Base {
  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    Base::tiePointers(layerData, dynRup, e_interoperability);

    auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*slipRateStrike)[size] = layerData.var(concreteLts->slipRateStrike);
    real(*slipRateDip)[size] = layerData.var(concreteLts->slipRateDip);
    real(*mu)[size] = layerData.var(concreteLts->mu);
    real(*stateVar)[size] = layerData.var(concreteLts->stateVariable);

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

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTOUT_RS_AGING_LAW_HPP
