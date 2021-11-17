#ifndef SEISSOL_DR_OUTOUT_LSW_BIMATERIAL_HPP
#define SEISSOL_DR_OUTOUT_LSW_BIMATERIAL_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public LinearSlipWeakening {
  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {
    LinearSlipWeakening::tiePointers(layerData, dynRup, e_interoperability);

    auto concreteLts =
        dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*regularisedStrength)[size] = layerData.var(concreteLts->regularisedStrength);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranStrength(
          ltsFace, meshFace, regularisedStrength);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTOUT_LSW_BIMATERIAL_HPP
