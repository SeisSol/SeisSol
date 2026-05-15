// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_SOLVER_SETTINGS_H_
#define SEISSOL_SRC_SOLVER_SETTINGS_H_

namespace seissol {
struct SimulationSettings {
  bool plasticity{false};
  bool integrate{false};

  SimulationSettings(bool plasticity, bool integrate)
      : plasticity(plasticity), integrate(integrate) {}
};
} // namespace seissol
#endif // SEISSOL_SRC_SOLVER_SETTINGS_H_
