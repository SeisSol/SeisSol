// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::geometry {
class MeshReader;
} // namespace seissol::geometry

namespace seissol::initializer {

struct ClusterLayout {
  std::vector<std::uint64_t> rates;
  double minimumTimestep;
  std::size_t globalClusterCount;

  ClusterLayout(const std::vector<std::uint64_t>& rates,
                double minimumTimestep,
                std::size_t globalClusterCount)
      : rates(rates), minimumTimestep(minimumTimestep), globalClusterCount(globalClusterCount) {}

  [[nodiscard]] double timestepRate(std::size_t id) const {
    return clusterRate(id) * minimumTimestep;
  }

  [[nodiscard]] std::uint64_t clusterRate(std::size_t id) const {
    std::uint64_t value = 1;
    for (std::size_t i = 0; i < id; ++i) {
      const auto rate = rates.size() > i ? rates[i] : rates.back();
      value *= rate;
    }
    return value;
  }

  static ClusterLayout fromMesh(const std::vector<std::uint64_t>& rates,
                                const geometry::MeshReader& mesh,
                                double wiggle,
                                bool infoprint);
};

} // namespace seissol::initializer
#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
