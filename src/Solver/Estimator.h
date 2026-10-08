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

#include <vector>

namespace seissol::solver {

/// The time the local kernel of `config` takes on this rank for a fixed amount of work.
auto miniSeisSol(ConfigId config) -> double;

/// How the work of a cell depends on its configuration: the cost of a cell of each configuration in
/// `configs` relative to a cell of `reference`, indexed by the id of the configuration; 1 for the
/// configurations built that are not in `configs`.
///
/// With `measure`, every rank times the time, local and neighbor kernels of each configuration with
/// the proxy, and all ranks take the median of the ratios over the ranks; collective then.
/// Otherwise, the ratios are those of the hardware FLOPs the kernels count, weighted with the size
/// of a real; that sees nothing of the machine.
auto configCostFactors(const std::vector<ConfigId>& configs, ConfigId reference, bool measure)
    -> std::vector<double>;

auto hostDeviceSwitch() -> int;

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_ESTIMATOR_H_
