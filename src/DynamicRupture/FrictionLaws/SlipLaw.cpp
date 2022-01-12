#include "SlipLaw.h"

namespace seissol::dr::friction_law {

real SlipLaw::updateStateVariable(int pointIndex,
                                  unsigned int face,
                                  real stateVarReference,
                                  real timeIncrement,
                                  real localSlipRate) {
  real localSl0 = sl0[face][pointIndex];
  real exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
  return localSl0 / localSlipRate * std::pow(localSlipRate * stateVarReference / localSl0, exp1);
}

} // namespace seissol::dr::friction_law
