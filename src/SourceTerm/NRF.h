// SPDX-FileCopyrightText: 2015 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_SOURCETERM_NRF_H_
#define SEISSOL_SRC_SOURCETERM_NRF_H_

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <vector>

namespace seissol::sourceterm {

// named like the compound type of NRF files that rconv writes with it
// NOLINTNEXTLINE(readability-identifier-naming)
struct Subfault_units {
  char* tinit;
  char* timestep;
  char* mu;
  char* area;
  char* tan1;
  char* tan2;
  char* normal;
};

using Subfault = struct Subfault {
  double tinit{};
  double timestep{};
  double mu{};
  double area{};
  Eigen::Vector3d tan1;
  Eigen::Vector3d tan2;
  Eigen::Vector3d normal;
};

using Offsets = std::array<unsigned, 3U>;

struct NRF {
  std::vector<Eigen::Vector3d> centres;
  std::vector<Subfault> subfaults;
  std::vector<Offsets> sroffsets;
  std::array<std::vector<double>, 3U> sliprates;
  [[nodiscard]] std::size_t size() const { return centres.size(); }
};
} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_NRF_H_
