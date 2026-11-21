// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_PLASTICITYPARAMETERS_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_PLASTICITYPARAMETERS_H_

#include <unordered_set>

namespace seissol::initializer::parameters {

struct PlasticityParameters {
  bool enabled{false};
  // the plasticity of every point of a cell, rather than one for the cell
  bool pointwise{true};
  // the mesh groups whose cells have none
  std::unordered_set<int> disabledGroups;
  // the relaxation time
  double tv{};
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_PLASTICITYPARAMETERS_H_
