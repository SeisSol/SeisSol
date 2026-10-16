// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_

#include "Initializer/Clustering/ClusterLadder.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace seissol::geometry {
class MeshReader;
} // namespace seissol::geometry

namespace seissol::initializer {

struct ClusterLayout {
  private:
  // declared first: `rates` is initialized from it
  ClusterLadder ladder_;

  public:
  /// Normalized ratios, i.e. exactly `globalClusterCount - 1` entries with the rate
  /// vector's abbreviations expanded.
  std::vector<std::uint64_t> rates;
  /// Base timestep; the wiggle factor is already applied.
  double minimumTimestep;
  std::size_t globalClusterCount;

  /// What fromMesh finds out about the clustering when it reports it, over all ranks.
  struct Report {
    /// the smallest cell timestep, before the wiggle factor is applied
    double minimumCellTimestep{0};
    double wiggleFactor{1};
    /// how many times fewer cell updates elementwise and clustered LTS take than global time
    /// stepping does
    double elementwiseSpeedup{1};
    double clusteredSpeedup{1};
    /// cells and dynamic rupture faces per cluster
    std::vector<std::uint64_t> cells;
    std::vector<std::uint64_t> dynamicRuptureFaces;
  };
  std::optional<Report> report;

  ClusterLayout(const std::vector<std::uint64_t>& rates,
                double minimumTimestep,
                std::size_t globalClusterCount)
      : ladder_(ClusterLadder::ofSize(rates, minimumTimestep, globalClusterCount)),
        rates(ladder_.ratios()), minimumTimestep(minimumTimestep),
        globalClusterCount(globalClusterCount) {}

  [[nodiscard]] double timestepRate(std::size_t id) const { return ladder_.timestep(id); }

  [[nodiscard]] std::uint64_t clusterRate(std::size_t id) const { return ladder_.updateFactor(id); }

  [[nodiscard]] const ClusterLadder& clusterLadder() const { return ladder_; }

  static ClusterLayout fromMesh(const std::vector<std::uint64_t>& rates,
                                const geometry::MeshReader& mesh,
                                double wiggle,
                                bool infoprint);
};

} // namespace seissol::initializer
#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
