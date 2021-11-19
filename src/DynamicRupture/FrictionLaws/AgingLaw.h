#ifndef SEISSOL_AGINGLAW_H
#define SEISSOL_AGINGLAW_H

#include "RateAndState.h"

namespace seissol::dr::friction_law {

/**
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class AgingLaw : public RateAndStateBase<AgingLaw> {
  public:
  using RateAndStateBase<AgingLaw>::RateAndStateBase;
  using RateAndStateBase<AgingLaw>::copyLtsTreeToLocal;

  protected:
  real calcStateVariableHook(real SV0, real tmp, real time_inc, real rs_sl0);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_AGINGLAW_H
