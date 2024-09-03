#ifndef SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
#define SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.h"

namespace seissol::dr::output {
class ImposedSlipRates : public ReceiverOutput {
  protected:
  real computeLocalStrength(LocalInfo& local) override { return 0.0; }

  void adjustRotatedUpdatedStress(std::array<real, 6>& rotatedUpdatedStress,
                                  const std::array<real, 6>& rotatedStress) override {
    // we plot the Stress from Godunov state, because we want
    // to see the traction change from the imposed slip distribution
    using namespace misc::quantity_indices;

    rotatedUpdatedStress[QuantityIndices::XY] = rotatedStress[QuantityIndices::XY];
    rotatedUpdatedStress[QuantityIndices::XZ] = rotatedStress[QuantityIndices::XZ];
  };
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_IMPOSED_RS_HPP
