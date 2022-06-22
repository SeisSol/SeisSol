#ifndef SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
#define SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class ImposedSlipRates : public ReceiverBasedOutput {
  protected:
  real computeLocalStrength() override { return 0.0; }

  void adjustRotatedUpdatedStress(std::array<real, 6>& rotatedUpdatedStress,
                                  std::array<real, 6>& rotatedStress) override {
    // we plot the Stress from Godunov state, because we want
    // to see the traction change from the imposed slip distribution
    rotatedUpdatedStress[3] = rotatedStress[3];
    rotatedUpdatedStress[5] = rotatedStress[5];
  };
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
