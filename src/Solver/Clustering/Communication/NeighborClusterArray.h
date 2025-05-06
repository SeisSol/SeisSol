// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTERARRAY_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTERARRAY_H_

#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <memory>
#include <vector>

#include "utils/logger.h"

namespace seissol::solver::clustering::communication {
template <typename ClusterT>
class NeighborClusterArray {
  public:
  NeighborClusterArray(std::vector<std::shared_ptr<ClusterT>>&& clusters, std::size_t rate)
      : clusters(std::move(clusters)), rate(rate) {}

  bool poll() {
    bool allTrue = true;
    std::size_t localCounter = counterStart;
    for (const auto& cluster : clusters) {
      if (!cluster->poll()) {
        allTrue = false;
      }
      if (localCounter % rate != 0) {
        break;
      }
      localCounter /= rate;
    }
    return allTrue;
  }

  void startFrom(parallel::runtime::StreamRuntime& runtime, bool sync) {
    ++counterStart;
    if (sync) {
      counterStart = 0;
    }

    std::size_t localCounter = counterStart;
    for (const auto& cluster : clusters) {
      cluster->startFrom(runtime);
      if (localCounter % rate != 0) {
        break;
      }
      localCounter /= rate;
    }
  }

  void stopTo(parallel::runtime::StreamRuntime& runtime, bool forceAll) {
    // ignored
    if (forceAll) {
      counterStart = 0;
    }
  }

  bool blocking() {
    return std::any_of(
        clusters.begin(), clusters.end(), [](const std::shared_ptr<ClusterT>& cluster) {
          return cluster->blocking();
        });
  }

  void dispose() {
    for (auto& cluster : clusters) {
      cluster->dispose();
    }
  }

  private:
  std::vector<std::shared_ptr<ClusterT>> clusters;
  std::size_t counterStart = 0;
  std::size_t rate = 0;
};
} // namespace seissol::solver::clustering::communication

#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTERARRAY_H_
