// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_IMPOSEDSLIPRATES_H_

#include "Common/Real.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::output {
class ImposedSlipRates : public ReceiverOutputImpl<ImposedSlipRates> {
  public:
  template <typename Cfg>
  Real<Cfg> computeLocalStrength(LocalInfo<Cfg>& /*local*/) {
    return 0.0;
  }

  /// The tractions that the imposed slip rate leaves, which no strength sets. Where the impedance
  /// couples the normal and the shear directions, they include the fault-normal traction.
  template <typename Cfg>
  void updateLocalTractions(LocalInfo<Cfg>& local,
                            Real<Cfg> /*strength*/,
                            Real<Cfg> /*strengthSlope*/) {
    applySlipRate(local,
                  stateAt<DynamicRupture::SlipRate1>(local, local.gpIndex),
                  stateAt<DynamicRupture::SlipRate2>(local, local.gpIndex));
  }

  template <typename Cfg>
  void adjustRotatedUpdatedStress(std::array<Real<Cfg>, 6>& rotatedUpdatedStress,
                                  const std::array<Real<Cfg>, 6>& rotatedStress) {
    // we plot the Stress from Godunov state, because we want
    // to see the traction change from the imposed slip distribution
    using namespace misc::quantity_indices;

    rotatedUpdatedStress[QuantityIndices::XY] = rotatedStress[QuantityIndices::XY];
    rotatedUpdatedStress[QuantityIndices::XZ] = rotatedStress[QuantityIndices::XZ];
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_IMPOSEDSLIPRATES_H_
