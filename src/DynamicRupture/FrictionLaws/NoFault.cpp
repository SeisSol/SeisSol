#include "NoFault.h"

namespace seissol::dr::friction_law {
template <typename Config>
void NoFault<Config>::updateFrictionAndSlip(
    FaultStresses<Config> const& faultStresses,
    TractionResults<Config>& tractionResults,
    std::array<RealT, misc::numPaddedPoints<Config>> const& stateVariableBuffer,
    std::array<RealT, misc::numPaddedPoints<Config>> const& strengthBuffer,
    unsigned ltsFace,
    unsigned timeIndex) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
    tractionResults.traction1[timeIndex][pointIndex] =
        faultStresses.traction1[timeIndex][pointIndex];
    tractionResults.traction2[timeIndex][pointIndex] =
        faultStresses.traction2[timeIndex][pointIndex];
  }
}
} // namespace seissol::dr::friction_law
