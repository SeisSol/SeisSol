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

  real computeLocalStrength() override {
    return -1.0 * local.mu * std::min(local.p + local.p0 - local.pf, static_cast<real>(0.0));
  }

  void adjustRotatedTractionAndStresses(std::array<real, 6>& rotatedTraction,
                                        std::array<real, 6>& rotatedLocalStress) override {
    rotatedTraction[3] = rotatedLocalStress[3];
    rotatedTraction[5] = rotatedLocalStress[5];
  };
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
