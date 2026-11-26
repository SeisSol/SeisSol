// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_

#include "Features/DynamicRupture/Solver/FrictionLaws/SourceTimeFunction.h"

namespace seissol::dr::friction_law::cpu {

// the source time functions are the same for the host and the device
template <typename Cfg>
using YoffeSTF = friction_law::YoffeSTF<Cfg>;
template <typename Cfg>
using GaussianSTF = friction_law::GaussianSTF<Cfg>;
template <typename Cfg>
using DeltaSTF = friction_law::DeltaSTF<Cfg>;

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SOURCETIMEFUNCTION_H_
