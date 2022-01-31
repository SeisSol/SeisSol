#include "SlipLaw.h"

namespace seissol::dr::friction_law {

double SlipLaw::updateStateVariable(int pointIndex,
                                  unsigned int face,
                                  double stateVarReference,
                                  double timeIncrement,
                                  double localSlipRate) {
  double localSl0 = sl0[face][pointIndex];
  double exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
  return localSl0 / localSlipRate * std::pow(localSlipRate * stateVarReference / localSl0, exp1);
}

} // namespace seissol::dr::friction_law
