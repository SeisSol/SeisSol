#include "NoFault.h"

namespace seissol::dr::friction_law {
void NoFault::updateFrictionAndSlip(FaultStresses& faultStresses,
                                    std::array<real, numPaddedPoints>& stateVariableBuffer,
                                    std::array<real, numPaddedPoints>& strengthBuffer,
                                    unsigned& ltsFace,
                                    unsigned& timeIndex) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    faultStresses.XYTractionResultGP[timeIndex][pointIndex] =
        faultStresses.XYStressGP[timeIndex][pointIndex];
    faultStresses.XZTractionResultGP[timeIndex][pointIndex] =
        faultStresses.XZStressGP[timeIndex][pointIndex];
  }
}
void NoFault::preHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace){};
void NoFault::postHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace){};
} // namespace seissol::dr::friction_law