#ifndef SEISSOL_FRICTIONLAWS_H
#define SEISSOL_FRICTIONLAWS_H

// collect all friction laws here
#include "AgingLaw.h"
#include "BaseFrictionLaw.h"
#include "FastVelocityWeakeningLaw.h"
#include "ImposedSlipRates.h"
#include "LinearSlipWeakening.h"
#include "NoFault.h"
#include "RateAndState.h"
#include "SourceTimeFunction.h"
#include "SlipLaw.h"
#include "SlowVelocityWeakeningLaw.h"

#ifdef ACL_DEVICE_OFFLOAD
#include "GpuImpl/LinearSlipWeakening.h"
#endif

#endif // SEISSOL_FRICTIONLAWS_H
