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
  /// Whether the material is sampled at the nodal points of each cell, so that
  /// a cell has to keep those samples rather than one value for itself.
  bool materialNodal{false};

  SimulationSettings(bool plasticity, bool integrate, bool materialNodal)
      : plasticity(plasticity), integrate(integrate), materialNodal(materialNodal) {}
};
} // namespace seissol
#endif // SEISSOL_SRC_SOLVER_SETTINGS_H_
