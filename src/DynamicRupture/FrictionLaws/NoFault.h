#ifndef SEISSOL_NOFAULT_H
#define SEISSOL_NOFAULT_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * No friction computation
 * input stress XYStressGP, XZStressGP equals output XYTractionResultGP, XZTractionResultGP
 */
class NoFault : public BaseFrictionLaw<NoFault> {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             std::array<real, numPaddedPoints>& stateVariableBuffer,
                             std::array<real, numPaddedPoints>& strengthBuffer,
                             unsigned& ltsFace,
                             unsigned& timeIndex);

  void preHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace);
  void postHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace);
};
} // namespace seissol::dr::friction_law
#endif // SEISSOL_NOFAULT_H
