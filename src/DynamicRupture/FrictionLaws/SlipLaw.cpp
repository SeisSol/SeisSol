#include "SlipLaw.h"

real seissol::dr::friction_law::SlipLaw::calcStateVariableHook(real SV0, real tmp, real time_inc, real RS_sl0) {
  return RS_sl0/tmp*seissol::dr::aux::power(tmp*SV0/RS_sl0, exp(-tmp*time_inc/RS_sl0));
}
