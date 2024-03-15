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
#include "SlipLaw.h"
#include "SlowVelocityWeakeningLaw.h"
#include "SourceTimeFunction.h"
#include "ThermalPressurization/NoTP.h"

#ifdef ACL_DEVICE
#include "GpuImpl/LinearSlipWeakening.h"
#include "GpuImpl/AgingLaw.h"
#include "GpuImpl/SlipLaw.h"
#include "GpuImpl/FastVelocityWeakeningLaw.h"
#include "GpuImpl/ThermalPressurization/NoTP.h"
#include "GpuImpl/ImposedSlipRates.h"
#endif

#endif // SEISSOL_FRICTIONLAWS_H
