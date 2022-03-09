#ifndef SEISSOL_DR_OUTPUT_RS_HPP
#define SEISSOL_DR_OUTPUT_RS_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class RateAndState : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* drDescr,
                   seissol::Interoperability& eInteroperability) override {
    Base::tiePointers(layerData, drDescr, eInteroperability);

    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(drDescr);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto numGaussPoints2d = init::QInterpolated::Stop[0];
    real(*slipRate1)[numGaussPoints2d] = layerData.var(concreteLts->slipRate1);
    real(*slipRate2)[numGaussPoints2d] = layerData.var(concreteLts->slipRate2);
    real(*mu)[numGaussPoints2d] = layerData.var(concreteLts->mu);
    real(*stateVar)[numGaussPoints2d] = layerData.var(concreteLts->stateVariable);
    real(*initialStressInFaultCS)[numGaussPoints2d][6] =
        layerData.var(concreteLts->initialStressInFaultCS);
    using VectorOfArraysT = std::vector<std::array<real, numGaussPoints2d>>;

    VectorOfArraysT initialStressXX(layerData.getNumberOfCells());
    VectorOfArraysT initialStressYY(layerData.getNumberOfCells());
    VectorOfArraysT initialStressZZ(layerData.getNumberOfCells());
    VectorOfArraysT initialStressXY(layerData.getNumberOfCells());
    VectorOfArraysT initialStressXZ(layerData.getNumberOfCells());
    VectorOfArraysT initialStressYZ(layerData.getNumberOfCells());

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      eInteroperability.copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRate1, slipRate2, mu);
      eInteroperability.copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVar);
      eInteroperability.copyFrictionOutputToFortranInitialStressInFaultCS(ltsFace,
                                                                          meshFace,
                                                                          initialStressInFaultCS,
                                                                          initialStressXX,
                                                                          initialStressYY,
                                                                          initialStressZZ,
                                                                          initialStressXY,
                                                                          initialStressYZ,
                                                                          initialStressXZ);
    }
  }

  void postCompute(seissol::initializers::DynamicRupture& drDescr) override {
    // do nothing
  }

  protected:
  real computeLocalStrength() override {
    return -1.0 * local.mu * std::min(local.p + local.p0 - local.pf, static_cast<real>(0.0));
  }
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_RS_HPP
