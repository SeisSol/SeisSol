// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

#include "NoFault.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "Kernels/Precision.h"
#include <array>

namespace seissol::dr::friction_law {
void NoFault::updateFrictionAndSlip(const FaultStresses& faultStresses,
                                    TractionResults& tractionResults,
                                    std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                                    std::array<real, misc::NumPaddedPoints>& strengthBuffer,
                                    unsigned ltsFace,
                                    unsigned timeIndex) {
  for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
    tractionResults.traction1[timeIndex][pointIndex] =
        faultStresses.traction1[timeIndex][pointIndex];
    tractionResults.traction2[timeIndex][pointIndex] =
        faultStresses.traction2[timeIndex][pointIndex];
  }
}
} // namespace seissol::dr::friction_law
