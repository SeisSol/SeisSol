#include "NoFault.h"

namespace seissol::dr::friction_law {
void NoFault::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                 seissol::initializers::DynamicRupture* dynRup,
                                 real fullUpdateTime) {}
void NoFault::updateFrictionAndSlip(FaultStresses& faultStresses,
                                    TractionResults& tractionResults,
                                    std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                                    std::array<real, misc::numPaddedPoints>& strengthBuffer,
                                    unsigned& ltsFace,
                                    unsigned& timeIndex) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    tractionResults.xyTraction[timeIndex][pointIndex] =
        faultStresses.xyStress[timeIndex][pointIndex];
    tractionResults.xzTraction[timeIndex][pointIndex] =
        faultStresses.xzStress[timeIndex][pointIndex];
  }
}
void NoFault::preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                      unsigned ltsFace){};
void NoFault::postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                       unsigned ltsFace){};
void NoFault::saveDynamicStressOutput(unsigned int ltsFace) {}
} // namespace seissol::dr::friction_law