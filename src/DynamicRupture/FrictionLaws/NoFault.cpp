#include "NoFault.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.hpp"
#include "Kernels/precision.hpp"
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
