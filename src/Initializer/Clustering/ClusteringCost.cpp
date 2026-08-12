// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Clustering/ClusteringCost.h"

#include "Initializer/Clustering/ClusterLadder.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <vector>

namespace seissol::initializer {

double computeLocalCostOfClustering(const std::vector<int>& clusterIds,
                                    const std::vector<std::uint64_t>& cellCosts,
                                    const std::vector<uint64_t>& rate,
                                    double wiggleFactor,
                                    double minimalTimestep) {
  assert(clusterIds.size() == cellCosts.size());

  double cost = 0.0;
  for (auto i = 0U; i < clusterIds.size(); ++i) {
    const auto cluster = clusterIds[i];
    const auto cellCost = cellCosts[i];
    const auto invUpdateFactor = ratepow(rate, 0, cluster);
    cost += static_cast<double>(cellCost) / static_cast<double>(invUpdateFactor);
  }

  const auto minDtWithWiggle = minimalTimestep * wiggleFactor;
  return cost / minDtWithWiggle;
}

double computeGlobalCostOfClustering(const std::vector<int>& clusterIds,
                                     const std::vector<std::uint64_t>& cellCosts,
                                     const std::vector<uint64_t>& rate,
                                     double wiggleFactor,
                                     double minimalTimestep,
                                     MPI_Comm comm) {
  double cost =
      computeLocalCostOfClustering(clusterIds, cellCosts, rate, wiggleFactor, minimalTimestep);
  MPI_Allreduce(MPI_IN_PLACE, &cost, 1, MPI_DOUBLE, MPI_SUM, comm);

  return cost;
}

std::vector<double>
    computeLocalCostsOfCappedClusterings(const std::vector<int>& clusterIds,
                                         const std::vector<std::uint64_t>& cellCosts,
                                         const std::vector<uint64_t>& rate,
                                         double wiggleFactor,
                                         double minimalTimestep,
                                         int maxClusterId) {
  assert(clusterIds.size() == cellCosts.size());
  assert(maxClusterId >= 0);

  const auto capCount = static_cast<std::size_t>(maxClusterId) + 1;

  // ratepow() rather than ClusterLadder: this is also reachable with hand-made clusterings
  // whose ids exceed what the rate vector could produce, which a ladder would reject.
  std::vector<double> updateFactors(capCount);
  for (std::size_t cluster = 0; cluster < capCount; ++cluster) {
    updateFactors[cluster] = static_cast<double>(ratepow(rate, 0, cluster));
  }

  std::vector<double> costs(capCount, 0.0);
  for (auto i = 0U; i < clusterIds.size(); ++i) {
    const auto cellCost = cellCosts[i];
    const auto cluster = clusterIds[i];
    for (std::size_t cap = 0; cap < capCount; ++cap) {
      const auto capped = std::min(static_cast<std::size_t>(cluster), cap);
      costs[cap] += static_cast<double>(cellCost) / updateFactors[capped];
    }
  }

  const auto minDtWithWiggle = minimalTimestep * wiggleFactor;
  for (auto& cost : costs) {
    cost /= minDtWithWiggle;
  }
  return costs;
}

std::vector<int> enforceMaxClusterId(const std::vector<int>& clusterIds, int maxClusterId) {
  auto newClusterIds = clusterIds;
  assert(maxClusterId >= 0);
  std::for_each(newClusterIds.begin(), newClusterIds.end(), [maxClusterId](auto& clusterId) {
    clusterId = std::min(clusterId, maxClusterId);
  });

  return newClusterIds;
}

// Merges clusters such that new cost is max oldCost * allowedPerformanceLossRatio
int computeMaxClusterIdAfterAutoMerge(const std::vector<int>& clusterIds,
                                      const std::vector<std::uint64_t>& cellCosts,
                                      const std::vector<uint64_t>& rate,
                                      double maximalAdmissibleCost,
                                      double wiggleFactor,
                                      double minimalTimestep) {
  int maxClusterId =
      clusterIds.empty() ? 0 : *std::max_element(clusterIds.begin(), clusterIds.end());
  MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, Mpi::mpi.comm());

  // All candidate costs in one pass, then one reduction -- previously this was a full
  // clustering copy plus an Allreduce for every candidate.
  auto costs = computeLocalCostsOfCappedClusterings(
      clusterIds, cellCosts, rate, wiggleFactor, minimalTimestep, maxClusterId);
  MPI_Allreduce(MPI_IN_PLACE,
                costs.data(),
                static_cast<int>(costs.size()),
                MPI_DOUBLE,
                MPI_SUM,
                Mpi::mpi.comm());

  // Iteratively merge clusters until we found the first number of clusters that has a cost that is
  // too high
  for (auto curMaxClusterId = maxClusterId; curMaxClusterId >= 0; --curMaxClusterId) {
    if (costs[curMaxClusterId] > maximalAdmissibleCost) {
      // This is the first number of clusters that resulted in an inadmissible cost
      // Hence, it was admissible in the previous iteration
      return std::min(maxClusterId, curMaxClusterId + 1);
    }
  }
  return 0;
}

std::uint64_t getCluster(double timestep,
                         double globalMinTimestep,
                         double ltsWiggleFactor,
                         const std::vector<std::uint64_t>& rate) {
  // Kept as a standalone loop rather than going through ClusterLadder: this is called for
  // single ad-hoc queries where building a ladder would allocate. It shares the ladder's
  // notion of where the rate vector terminates, and the unit tests pin the two together.
  if (rate.empty()) {
    return 0;
  }

  const auto clusterLimit = ClusterLadder::intrinsicClusterCount(rate);

  double upper = ltsWiggleFactor * rate[0] * globalMinTimestep;

  std::uint64_t cluster = 0;
  while (upper <= timestep) {
    if (clusterLimit != ClusterLadder::Unbounded && cluster + 1 >= clusterLimit) {
      break;
    }
    const auto currentRate = rate.size() > (cluster + 1) ? rate[cluster + 1] : rate.back();
    upper *= currentRate;
    ++cluster;
  }
  return cluster;
}

std::uint64_t ratepow(const std::vector<std::uint64_t>& rate, std::uint64_t a, std::uint64_t b) {
  std::uint64_t factor = 1;
  for (std::uint64_t i = a; i < b; ++i) {
    factor *= i < rate.size() ? rate[i] : rate.back();
  }
  return factor;
}

} // namespace seissol::initializer
