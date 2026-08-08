// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERHISTOGRAM_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERHISTOGRAM_H_

#include "Initializer/TimeStepping/ClusterLadder.h"

#include <cstddef>
#include <limits>
#include <mpi.h>
#include <vector>

namespace seissol::initializer {

/// Cell cost summed per time cluster.
///
/// Every clustering cost the LTS search needs is a function of this histogram and a
/// ClusterLadder alone: once it is built, neither the mesh nor the per-cell arrays are
/// touched again, and no further communication is required. That is what makes candidate
/// evaluation local and O(clusterCount) instead of O(cells) plus a reduction.
///
/// # Not a drop-in replacement for the per-cell cost
///
/// `cost()` sums `weight(k) / updateFactor(k)`, whereas `computeLocalCostOfClustering()`
/// sums `cellCost_i / updateFactor(clusterId_i)` in cell order. For rate vectors built from
/// powers of two every term is a dyadic rational and the two agree bit for bit; for other
/// rates they drift by a few ulp (measured: ~1e-13 relative at 3e6 cells with rate 3).
///
/// The legacy search therefore keeps using the per-cell form, and this type is the
/// substrate for the new one rather than a replacement underneath the old one.
class ClusterHistogram {
  public:
  static constexpr std::size_t NoCap = std::numeric_limits<std::size_t>::max();

  ClusterHistogram() = default;
  explicit ClusterHistogram(std::vector<double> weights);

  /// Accumulates `cellCosts` into `clusterCount` bins indexed by `clusterIds`.
  static ClusterHistogram fromClustering(const std::vector<int>& clusterIds,
                                         const std::vector<int>& cellCosts,
                                         std::size_t clusterCount);

  /// Sums the histogram across all ranks in place. The only collective the search needs.
  void reduce(MPI_Comm comm);

  [[nodiscard]] std::size_t clusterCount() const { return weights_.size(); }
  [[nodiscard]] const std::vector<double>& weights() const { return weights_; }
  [[nodiscard]] double weight(std::size_t cluster) const { return weights_[cluster]; }
  [[nodiscard]] double totalWeight() const;

  /// Cost per unit of simulated time. Clusters above `maxCluster` are folded into it, which
  /// is the histogram equivalent of `enforceMaxClusterId()`.
  [[nodiscard]] double cost(const ClusterLadder& ladder, std::size_t maxCluster = NoCap) const;

  private:
  std::vector<double> weights_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERHISTOGRAM_H_
