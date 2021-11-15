#ifndef SEISSOL_SLIPLAW_H
#define SEISSOL_SLIPLAW_H

#include "AgingLaw.h"

namespace seissol::dr::friction_law {
/**
 * This class was not tested and compared to the Fortran FL4. Since FL4 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class SlipLaw : public AgingLaw {
  public:
  using AgingLaw::AgingLaw;

  virtual real calcStateVariableHook(real SV0, real tmp, real time_inc, real RS_sl0) override;
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_SLIPLAW_H
