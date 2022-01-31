#include "AgingLaw.h"
namespace seissol::dr::friction_law {
double AgingLaw::updateStateVariable(int pointIndex,
                                   unsigned int face,
                                   double stateVarReference,
                                   double timeIncrement,
                                   double localSlipRate) {
  double localSl0 = sl0[face][pointIndex];
  double exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
  return stateVarReference * exp1 + localSl0 / localSlipRate * (1.0 - exp1);
}

} // namespace seissol::dr::friction_law