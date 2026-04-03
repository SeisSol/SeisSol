// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "NoFault.h"

#include "Common/Executor.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "Kernels/Precision.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {
void NoFault::updateFrictionAndSlip(
    const FaultStresses<Executor::Host>& faultStresses,
    TractionResults<Executor::Host>& tractionResults,
    std::array<real, misc::NumPaddedPoints>& /*stateVariableBuffer*/,
    std::array<real, misc::NumPaddedPoints>& /*strengthBuffer*/,
    std::size_t /*ltsFace*/,
    uint32_t timeIndex) {
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
    tractionResults.traction1[timeIndex][pointIndex] =
        faultStresses.traction1[timeIndex][pointIndex];
    tractionResults.traction2[timeIndex][pointIndex] =
        faultStresses.traction2[timeIndex][pointIndex];
  }
}
} // namespace seissol::dr::friction_law::cpu
