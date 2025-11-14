// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: David Schneller

#ifndef SEISSOL_SRC_SOURCETERM_FSRMREADER_H_
#define SEISSOL_SRC_SOURCETERM_FSRMREADER_H_

#include "Kernels/Precision.h"

#include <Eigen/Dense>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::sourceterm {

// TODO: when refactoring, replace raw array types
struct FSRMSource {
  double momentTensor[3][3]{};
  double solidVelocityComponent[3]{};
  double pressureComponent{};
  double fluidVelocityComponent[3]{};
  size_t numberOfSources{};
  std::vector<Eigen::Vector3d> centers;
  std::vector<double> strikes;
  std::vector<double> dips;
  std::vector<double> rakes;
  std::vector<double> onsets;
  std::vector<double> areas;
  double timestep{};
  size_t numberOfSamples{};
  std::vector<std::vector<double>> timeHistories;

  void read(const std::string& filename);
};
} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_FSRMREADER_H_
