// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_

#include "Initializer/TimeStepping/ClusterLadder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::geometry {
class MeshReader;
} // namespace seissol::geometry

namespace seissol::initializer {

struct ClusterLayout {
  private:
  // declared first: `rates` is initialized from it
  ClusterLadder ladder;

  public:
  /// Normalized ratios, i.e. exactly `globalClusterCount - 1` entries with the rate
  /// vector's abbreviations expanded.
  std::vector<std::uint64_t> rates;
  /// Base timestep; the wiggle factor is already applied.
  double minimumTimestep;
  std::size_t globalClusterCount;

  ClusterLayout(const std::vector<std::uint64_t>& rates,
                double minimumTimestep,
                std::size_t globalClusterCount)
      : ladder(ClusterLadder::ofSize(rates, minimumTimestep, globalClusterCount)),
        rates(ladder.ratios()), minimumTimestep(minimumTimestep),
        globalClusterCount(globalClusterCount) {}

  [[nodiscard]] double timestepRate(std::size_t id) const { return ladder.timestep(id); }

  [[nodiscard]] std::uint64_t clusterRate(std::size_t id) const { return ladder.updateFactor(id); }

  [[nodiscard]] const ClusterLadder& clusterLadder() const { return ladder; }

  static ClusterLayout fromMesh(const std::vector<std::uint64_t>& rates,
                                const geometry::MeshReader& mesh,
                                double wiggle,
                                bool infoprint);
};

} // namespace seissol::initializer
#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERLAYOUT_H_
