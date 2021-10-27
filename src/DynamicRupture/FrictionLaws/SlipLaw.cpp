#include "SlipLaw.h"

namespace seissol::dr::friction_law {

real SlipLaw::calcStateVariableHook(real SV0, real tmp, real time_inc, real RS_sl0) {
  return RS_sl0 / tmp * std::pow(tmp * SV0 / RS_sl0, exp(-tmp * time_inc / RS_sl0));
}
} // namespace seissol::dr::friction_law
