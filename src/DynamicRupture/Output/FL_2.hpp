#ifndef SEISSOL_DROUTOUT_FL_2_HPP
#define SEISSOL_DROUTOUT_FL_2_HPP

#include "DynamicRupture/Output/Base.hpp"

namespace seissol::dr::output {
class FL_2 : public Base {
  public:
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* dynRup,
                           seissol::Interoperability& e_interoperability) override {

    Base::tiePointers(layerData, dynRup, e_interoperability);

    seissol::initializers::DR_FL_2* ConcreteLts =
        dynamic_cast<seissol::initializers::DR_FL_2*>(dynRup);

    DRFaceInformation* faceInformation = layerData.var(ConcreteLts->faceInformation);
    real* averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
    real(*dynStress_time)[init::QInterpolated::Stop[0]] =
        layerData.var(ConcreteLts->dynStress_time);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averaged_Slip, dynStress_time);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture& DynRup) override {
    std::cout << "output vars for FL_2\n";
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DROUTOUT_FL_2_HPP
