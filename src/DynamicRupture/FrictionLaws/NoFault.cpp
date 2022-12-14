#include "NoFault.h"

namespace seissol::dr::friction_law {
void NoFault::updateFrictionAndSlip(
    FaultStresses const& faultStresses,
    TractionResults& tractionResults,
    std::array<real, misc::numPaddedPoints> const& stateVariableBuffer,
    std::array<real, misc::numPaddedPoints> const& strengthBuffer,
    unsigned ltsFace,
    unsigned timeIndex) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    tractionResults.traction1[timeIndex][pointIndex] =
        faultStresses.traction1[timeIndex][pointIndex];
    tractionResults.traction2[timeIndex][pointIndex] =
        faultStresses.traction2[timeIndex][pointIndex];
  }
}
} // namespace seissol::dr::friction_law
