// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: David Schneller

#ifndef SEISSOL_SRC_SOLVER_ESTIMATOR_H_
#define SEISSOL_SRC_SOLVER_ESTIMATOR_H_

#include "Common/ConfigRegistry.h"

namespace seissol::solver {

/// The time the local kernel of `config` takes on this rank for a fixed amount of work.
auto miniSeisSol(ConfigId config) -> double;

auto hostDeviceSwitch() -> int;

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_ESTIMATOR_H_
