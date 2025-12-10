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
template <typename Cfg>
void NoFault<Cfg>::updateFrictionAndSlip(
    const FaultStresses<Cfg, Executor::Host>& faultStresses,
    TractionResults<Cfg, Executor::Host>& tractionResults,
    std::array<real, misc::NumPaddedPoints<Cfg>>& /*stateVariableBuffer*/,
    std::array<real, misc::NumPaddedPoints<Cfg>>& /*strengthBuffer*/,
    std::size_t /*ltsFace*/,
    uint32_t timeIndex) {
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; pointIndex++) {
    tractionResults.traction1[timeIndex][pointIndex] =
        faultStresses.traction1[timeIndex][pointIndex];
    tractionResults.traction2[timeIndex][pointIndex] =
        faultStresses.traction2[timeIndex][pointIndex];
  }
}

#define SEISSOL_CONFIGITER(cfg) template class NoFault<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::dr::friction_law::cpu
