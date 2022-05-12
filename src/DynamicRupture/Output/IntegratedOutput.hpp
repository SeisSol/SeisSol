#ifndef SEISSOL_INTEGRATEDOUTPUT_HPP
#define SEISSOL_INTEGRATEDOUTPUT_HPP

#include "DynamicRupture/Output/ParametersInitializer.hpp"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/LTS.h"
#include "Initializer/DynamicRupture.h"
#include "Solver/Interoperability.h"

namespace seissol::dr::output {
class IntegratedOutput {
  public:
  void setLtsData(seissol::initializers::DynamicRupture* userDrDescr, size_t numFaultElements);
  void setFaceToLtsMap(FaceToLtsMapT* map) { faceToLtsMap = map; }

  long double getSeismicMoment(IntegratedOutputData& outputData);
  double getSeismicMomentRate(IntegratedOutputData& outputData);

  private:
  seissol::initializers::DynamicRupture* drDescr{nullptr};
  FaceToLtsMapT* faceToLtsMap{nullptr};
  size_t numFaultElements{0};
};
} // namespace seissol::dr::output

#endif // SEISSOL_INTEGRATEDOUTPUT_HPP
