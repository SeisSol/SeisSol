// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Clustering/ClusterHistogram.h"

#include "Initializer/Clustering/ClusterLadder.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <utility>
#include <vector>

namespace seissol::initializer {

ClusterHistogram::ClusterHistogram(std::vector<double> weights) : weights_(std::move(weights)) {}

ClusterHistogram ClusterHistogram::fromClustering(const std::vector<int>& clusterIds,
                                                  const std::vector<std::uint64_t>& cellCosts,
                                                  std::size_t clusterCount) {
  assert(clusterIds.size() == cellCosts.size());
  std::vector<double> weights(clusterCount, 0.0);
  for (std::size_t cell = 0; cell < clusterIds.size(); ++cell) {
    assert(clusterIds[cell] >= 0);
    const auto cluster = static_cast<std::size_t>(clusterIds[cell]);
    assert(cluster < clusterCount);
    weights[cluster] += static_cast<double>(cellCosts[cell]);
  }
  return ClusterHistogram(std::move(weights));
}

void ClusterHistogram::reduce(MPI_Comm comm) {
  if (weights_.empty()) {
    return;
  }
  MPI_Allreduce(
      MPI_IN_PLACE, weights_.data(), static_cast<int>(weights_.size()), MPI_DOUBLE, MPI_SUM, comm);
}

double ClusterHistogram::totalWeight() const {
  double total = 0.0;
  for (const auto weight : weights_) {
    total += weight;
  }
  return total;
}

double ClusterHistogram::cost(const ClusterLadder& ladder, std::size_t maxCluster) const {
  if (weights_.empty()) {
    return 0.0;
  }
  const auto cap = std::min(maxCluster, weights_.size() - 1);
  assert(cap < ladder.clusterCount());

  double cost = 0.0;
  for (std::size_t cluster = 0; cluster < cap; ++cluster) {
    cost += weights_[cluster] / static_cast<double>(ladder.updateFactor(cluster));
  }
  // everything above the cap is demoted onto it, as enforceMaxClusterId() would do
  double folded = 0.0;
  for (std::size_t cluster = cap; cluster < weights_.size(); ++cluster) {
    folded += weights_[cluster];
  }
  cost += folded / static_cast<double>(ladder.updateFactor(cap));

  return cost / ladder.baseTimestep();
}

} // namespace seissol::initializer
