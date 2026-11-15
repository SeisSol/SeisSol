// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_SOLVER_SETTINGS_H_
#define SEISSOL_SRC_SOLVER_SETTINGS_H_

#include <cstddef>

namespace seissol {
struct SimulationSettings {
  bool plasticity{false};
  bool integrate{false};

  // the number of simulations of the run; unlike the number of fused simulations, it is not fixed
  // by the build
  std::size_t simulationCount{1};

  // the number of simulations fused into a cell (the same for all configurations of a run)
  std::size_t fusedCount{1};

  SimulationSettings(bool plasticity,
                     bool integrate,
                     std::size_t simulationCount,
                     std::size_t fusedCount)
      : plasticity(plasticity), integrate(integrate), simulationCount(simulationCount),
        fusedCount(fusedCount) {}

  // the copies of every cell that hold the simulations: the copy c holds the simulations
  // c * fusedCount to (c + 1) * fusedCount - 1. The fused simulations of the last copy past
  // simulationCount only fill it up.
  [[nodiscard]] std::size_t copyCount() const {
    return (simulationCount + fusedCount - 1) / fusedCount;
  }
};
} // namespace seissol
#endif // SEISSOL_SRC_SOLVER_SETTINGS_H_
