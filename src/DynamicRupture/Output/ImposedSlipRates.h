// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_IMPOSEDSLIPRATES_H_

#include "DynamicRupture/Output/ReceiverBasedOutput.h"

namespace seissol::dr::output {
class ImposedSlipRates : public ReceiverOutput {
  protected:
  std::array<real, seissol::multisim::NumSimulations> computeLocalStrength(LocalInfo& local) override { return std::array<real, seissol::multisim::NumSimulations>{0.0};}

  void adjustRotatedUpdatedStress(std::array<std::array<real, 6>, seissol::multisim::NumSimulations>& rotatedUpdatedStress,
                                  const std::array<std::array<real, 6>, seissol::multisim::NumSimulations>& rotatedStress) override {
    // we plot the Stress from Godunov state, because we want
    // to see the traction change from the imposed slip distribution
    using namespace misc::quantity_indices;
    for(unsigned int sim=0; sim < seissol::multisim::NumSimulations; ++sim) {
    rotatedUpdatedStress[sim][QuantityIndices::XY] = rotatedStress[sim][QuantityIndices::XY];
    rotatedUpdatedStress[sim][QuantityIndices::XZ] = rotatedStress[sim][QuantityIndices::XZ];
    }
  };
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_IMPOSEDSLIPRATES_H_
