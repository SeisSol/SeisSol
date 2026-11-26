// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_

#include "Features/DynamicRupture/Solver/FrictionLaws/ThermalPressurization/NoTP.h"

namespace seissol::dr::friction_law::cpu {

/// the absence of thermal pressurization of the rate-and-state laws, the same on the host and the
/// device
template <typename Cfg>
using NoTP = friction_law::NoTP<Cfg>;

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_
