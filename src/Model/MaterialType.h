// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_MATERIALTYPE_H_
#define SEISSOL_SRC_MODEL_MATERIALTYPE_H_

namespace seissol::model {
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

#endif // SEISSOL_SRC_MODEL_MATERIALTYPE_H_
