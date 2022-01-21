#include "AgingLaw.h"
namespace seissol::dr::friction_law {
real AgingLaw::updateStateVariable(int pointIndex,
                                   unsigned int face,
                                   real stateVarReference,
                                   real timeIncrement,
                                   real localSlipRate) {
  real localSl0 = sl0[face][pointIndex];
  real exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
  return stateVarReference * exp1 + localSl0 / localSlipRate * (1.0 - exp1);
}

} // namespace seissol::dr::friction_law