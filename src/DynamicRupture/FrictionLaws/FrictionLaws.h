// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONLAWS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONLAWS_H_

// IWYU pragma: begin_exports

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
#include "GpuImpl/AgingLaw.h"
#include "GpuImpl/FastVelocityWeakeningLaw.h"
#include "GpuImpl/LinearSlipWeakening.h"
#include "GpuImpl/NoFault.h"
#include "GpuImpl/SlipLaw.h"
#include "GpuImpl/ThermalPressurization/NoTP.h"
#endif

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONLAWS_H_
