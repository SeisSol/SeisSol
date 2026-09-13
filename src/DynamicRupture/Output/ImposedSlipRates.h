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
  real computeLocalStrength(LocalInfo& /*local*/) override { return 0.0; }

  void adjustRotatedUpdatedStress(std::array<real, 6>& rotatedUpdatedStress,
                                  const std::array<real, 6>& rotatedStress) override {
    // we plot the Stress from Godunov state, because we want
    // to see the traction change from the imposed slip distribution
    rotatedUpdatedStress[misc::voigt::XY] = rotatedStress[misc::voigt::XY];
    rotatedUpdatedStress[misc::voigt::XZ] = rotatedStress[misc::voigt::XZ];
  };
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_IMPOSEDSLIPRATES_H_
