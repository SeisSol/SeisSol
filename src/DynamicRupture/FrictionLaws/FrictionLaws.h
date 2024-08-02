#ifndef SEISSOL_FRICTIONLAWS_H
#define SEISSOL_FRICTIONLAWS_H

// IWYU pragma: begin_exports

// collect all friction laws here
#include "AdjointRSF.h"
#include "AdjointSlowVelWeakening.h"
#include "AdjointSlip.h"
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
#include "GpuImpl/AgingLaw.h"
#include "GpuImpl/FastVelocityWeakeningLaw.h"
#include "GpuImpl/LinearSlipWeakening.h"
#include "GpuImpl/SlipLaw.h"
#include "GpuImpl/ThermalPressurization/NoTP.h"
#endif

// IWYU pragma: end_exports

#endif // SEISSOL_FRICTIONLAWS_H
