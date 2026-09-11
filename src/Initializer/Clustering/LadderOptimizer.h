// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_LADDEROPTIMIZER_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_LADDEROPTIMIZER_H_

#include "Initializer/Clustering/ClusterCostModel.h"
#include "Initializer/Clustering/TimestepHistogram.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::initializer {

struct LadderConstraints {
  /// Coarsest update factor any cluster may have. Bounded from above by the largest cell
  /// timestep, and in practice also by the synchronization interval.
  std::size_t maxIndex{1};
  /// Upper bound on the number of clusters.
  std::size_t maxClusterCount{64};
  /// Upper bound on a single ratio between adjacent clusters; 0 means unbounded.
  std::uint64_t maxRatio{0};
};

struct LadderCandidate {
  /// Normalized ratios, ready to hand to ClusterLadder. Never contains a 1.
  std::vector<std::uint64_t> ratios;
  double cost{0.0};
};

/// Finds the cheapest ladder for a given distribution of cell timesteps.
///
/// The update factors of a ladder form a divisibility chain `1 | n_1 | n_2 | ...`, because a
/// cluster may only advance when all finer clusters have caught up. Minimizing over all such
/// chains is a shortest path on the divisor lattice of `[1, maxIndex]`: the state is the
/// update factor of the current rung, and an edge to a proper multiple pays for the cells
/// that the current rung will carry.
///
/// Cost is `O(maxClusterCount * maxIndex * log(maxIndex))` and needs neither the mesh nor
/// communication once the histogram exists.
///
/// Ties are resolved towards fewer clusters. This matters: under the pure update-count
/// objective an empty rung is free, so without the tie-break the result would be padded with
/// clusters that carry nothing.
LadderCandidate optimalLadder(const TimestepHistogram& histogram,
                              const ClusterCostModel& costModel,
                              double baseTimestep,
                              const LadderConstraints& constraints);

/// Cost of one explicitly given ladder, for comparing against what the optimizer picked.
double ladderCost(const TimestepHistogram& histogram,
                  const ClusterCostModel& costModel,
                  double baseTimestep,
                  const std::vector<std::uint64_t>& ratios);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_LADDEROPTIMIZER_H_
