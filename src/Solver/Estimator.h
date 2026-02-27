// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: David Schneller

#ifndef SEISSOL_SRC_SOLVER_ESTIMATOR_H_
#define SEISSOL_SRC_SOLVER_ESTIMATOR_H_

namespace seissol::solver {

auto miniSeisSol() -> double;

auto hostDeviceSwitch() -> int;

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_ESTIMATOR_H_
