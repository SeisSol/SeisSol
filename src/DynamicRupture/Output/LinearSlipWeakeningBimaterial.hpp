#ifndef SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP
#define SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public LinearSlipWeakening {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* drDescr,
                   seissol::Interoperability& eInteroperability) override {
    LinearSlipWeakening::tiePointers(layerData, drDescr, eInteroperability);

    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(drDescr);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*regularisedStrength)[size] = layerData.var(concreteLts->regularisedStrength);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      eInteroperability.copyFrictionOutputToFortranStrength(ltsFace, meshFace, regularisedStrength);
    }
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_LSW_BIMATERIAL_HPP
