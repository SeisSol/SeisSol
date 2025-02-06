// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: David Schneller

#ifndef SEISSOL_SRC_SOURCETERM_FSRMREADER_H_
#define SEISSOL_SRC_SOURCETERM_FSRMREADER_H_

#include <Eigen/Dense>
#include <cstddef>
#include <string>
#include <vector>

#include "Kernels/Precision.h"

namespace seissol::sourceterm {

// TODO: when refactoring, replace raw array types
struct FSRMSource {
  real momentTensor[3][3]{};
  real solidVelocityComponent[3]{};
  real pressureComponent{};
  real fluidVelocityComponent[3]{};
  size_t numberOfSources{};
  std::vector<Eigen::Vector3d> centers;
  std::vector<real> strikes;
  std::vector<real> dips;
  std::vector<real> rakes;
  std::vector<real> onsets;
  std::vector<real> areas;
  real timestep{};
  size_t numberOfSamples{};
  std::vector<std::vector<real>> timeHistories;

  void read(const std::string& filename);
};
} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_FSRMREADER_H_
