// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONLAWS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONLAWS_H_

// IWYU pragma: begin_exports

// collect all friction laws here
#include "CpuImpl/AgingLaw.h"
#include "CpuImpl/BaseFrictionLaw.h"
#include "CpuImpl/FastVelocityWeakeningLaw.h"
#include "CpuImpl/ImposedSlipRates.h"
#include "CpuImpl/LinearSlipWeakening.h"
#include "CpuImpl/NoFault.h"
#include "CpuImpl/RateAndState.h"
#include "CpuImpl/SevereVelocityWeakeningLaw.h"
#include "CpuImpl/SlipLaw.h"
#include "CpuImpl/SlowVelocityWeakeningLaw.h"
#include "CpuImpl/SourceTimeFunction.h"
#include "CpuImpl/ThermalPressurization/NoTP.h"
#include "CpuImpl/ThermalPressurization/ThermalPressurization.h"

#ifdef ACL_DEVICE
#include "GpuImpl/AgingLaw.h"
#include "GpuImpl/FastVelocityWeakeningLaw.h"
#include "GpuImpl/ImposedSlipRates.h"
#include "GpuImpl/LinearSlipWeakening.h"
#include "GpuImpl/NoFault.h"
#include "GpuImpl/SevereVelocityWeakeningLaw.h"
#include "GpuImpl/SlipLaw.h"
#include "GpuImpl/SourceTimeFunction.h"
#include "GpuImpl/ThermalPressurization/NoTP.h"
#include "GpuImpl/ThermalPressurization/ThermalPressurization.h"
#endif

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONLAWS_H_
