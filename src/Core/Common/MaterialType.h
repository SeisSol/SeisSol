// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_COMMON_MATERIALTYPE_H_
#define SEISSOL_SRC_CORE_COMMON_MATERIALTYPE_H_

namespace seissol::model {

/// The kinds of material a configuration computes in; part of the vocabulary of the
/// configurations, as the solvers are (SolverType). What a kind of material is, the wave
/// propagation says (model::MaterialOf).
enum class MaterialType {
  Solid,
  Acoustic,
  Elastic,
  Viscoelastic,
  Viscoacoustic,
  Anisotropic,
  Poroelastic
};
} // namespace seissol::model

#endif // SEISSOL_SRC_CORE_COMMON_MATERIALTYPE_H_
