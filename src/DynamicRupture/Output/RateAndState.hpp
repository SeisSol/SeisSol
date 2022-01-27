#ifndef SEISSOL_DR_OUTPUT_RS_HPP
#define SEISSOL_DR_OUTPUT_RS_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class RateAndState : public Base {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* dynRup,
                   seissol::Interoperability& eInteroperability) override {
    Base::tiePointers(layerData, dynRup, eInteroperability);

    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto numGaussPoints2d = init::QInterpolated::Stop[0];
    real(*slipRate1)[numGaussPoints2d] = layerData.var(concreteLts->slipRate1);
    real(*slipRate2)[numGaussPoints2d] = layerData.var(concreteLts->slipRate2);
    real(*mu)[numGaussPoints2d] = layerData.var(concreteLts->mu);
    real(*stateVar)[numGaussPoints2d] = layerData.var(concreteLts->stateVariable);
    real(*initialStressInFaultCS)[numGaussPoints2d][6] =
        layerData.var(concreteLts->initialStressInFaultCS);
    using VectorOfArrays = std::vector<std::array<real, numGaussPoints2d>>;

    VectorOfArrays initialStressXX(layerData.getNumberOfCells());
    VectorOfArrays initialStressYY(layerData.getNumberOfCells());
    VectorOfArrays initialStressZZ(layerData.getNumberOfCells());
    VectorOfArrays initialStressXY(layerData.getNumberOfCells());
    VectorOfArrays initialStressXZ(layerData.getNumberOfCells());
    VectorOfArrays initialStressYZ(layerData.getNumberOfCells());

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

  void postCompute(seissol::initializers::DynamicRupture& dynRup) override {
    // do nothing
  }
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_RS_HPP
