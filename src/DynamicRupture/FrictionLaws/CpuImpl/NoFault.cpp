// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "NoFault.h"

#include "Common/Executor.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {
template <typename Cfg>
void NoFault<Cfg>::updateFrictionAndSlip(
    const FaultStresses<Cfg, Executor::Host>& faultStresses,
    const FaultStresses<Cfg, Executor::Host>& /*initialStress*/,
    TractionResults<Cfg, Executor::Host>& tractionResults,
    std::array<real, misc::NumPaddedPoints<Cfg>>& /*stateVariableBuffer*/,
    std::array<real, misc::NumPaddedPoints<Cfg>>& /*strengthBuffer*/,
    std::size_t /*ltsFace*/,
    uint32_t /*timeIndex*/) {
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; pointIndex++) {
    tractionResults.traction1[pointIndex] = faultStresses.traction1[pointIndex];
    tractionResults.traction2[pointIndex] = faultStresses.traction2[pointIndex];
  }
}
#define SEISSOL_INSTANTIATE(Cfg) template class NoFault<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::dr::friction_law::cpu
