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
  /// Values per cell of the state of the derived outputs of the wave field.
  std::size_t derivedState{0};

  SimulationSettings(bool plasticity, bool integrate)
      : plasticity(plasticity), integrate(integrate) {}
};
} // namespace seissol
#endif // SEISSOL_SRC_SOLVER_SETTINGS_H_
