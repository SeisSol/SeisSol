// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_

#include <vector>

namespace seissol::initializer {

struct ClusterLayout {
  std::vector<std::size_t> rates;
  double minimumTimestep;
  std::size_t globalClusterCount;

  ClusterLayout(const std::vector<std::size_t>& rates,
                double minimumTimestep,
                std::size_t globalClusterCount)
      : rates(rates), minimumTimestep(minimumTimestep), globalClusterCount(globalClusterCount) {}

  [[nodiscard]] double timestepRate(std::size_t id) const {
    return clusterRate(id) * minimumTimestep;
  }

  [[nodiscard]] long clusterRate(std::size_t id) const {
    std::size_t value = 1;
    for (std::size_t i = 0; i < id; ++i) {
      const auto rate = rates.size() > i ? rates[i] : rates.back();
      value *= rate;
    }
    return value;
  }
};

} // namespace seissol::initializer
#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
