#ifndef SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
#define SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class ImposedSlipRates : public ReceiverBasedOutput {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* drDescr,
                   seissol::Interoperability& eInteroperability) override {
    ReceiverBasedOutput::tiePointers(layerData, drDescr, eInteroperability);
  }

  real computeLocalStrength() override { return 0.0; }

  void adjustRotatedUpdatedStress(std::array<real, 6>& rotatedUpdatedStress,
                                  std::array<real, 6>& rotatedStress) override {
    rotatedUpdatedStress[3] = rotatedStress[3];
    rotatedUpdatedStress[5] = rotatedStress[5];
  };
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
