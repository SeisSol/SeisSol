#ifndef SEISSOL_NOFAULT_H
#define SEISSOL_NOFAULT_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * No friction computation input stress equals output
 */
class NoFault : public BaseFrictionLaw<NoFault> {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  static void updateFrictionAndSlip(const FaultStresses& faultStresses,
                                    TractionResults& tractionResults,
                                    std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                                    std::array<real, misc::NumPaddedPoints>& strengthBuffer,
                                    unsigned ltsFace,
                                    unsigned timeIndex);

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {};
  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {};
  void saveDynamicStressOutput(unsigned int ltsFace) {};
};
} // namespace seissol::dr::friction_law
#endif // SEISSOL_NOFAULT_H
