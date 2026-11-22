// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_FACECLUSTERFACTORY_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_FACECLUSTERFACTORY_H_

#include "Core/Common/Executor.h"
#include "Core/Monitoring/ActorStateStatistics.h"
#include "Core/Monitoring/DeviceTimeline.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/TimeStepping/Compute/FaceCluster.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"

#include <cstddef>
#include <memory>

namespace seissol::solver {

/// Where a cluster reports what it does.
struct ClusterMonitoring {
  LoopStatistics* loopStatistics{nullptr};
  DeviceTimeline* deviceTimeline{nullptr};
  ActorStateStatistics* actorStateStatistics{nullptr};
};

/**
 * Makes the clusters of a feature whose work lies on faces between cells, such as the dynamic
 * rupture: at most one per layer of the cells, as the layers of its faces have the colors of the
 * cells. The faces of an interior layer lie between interior or copy cells, those of a copy layer
 * between copy and ghost cells; the time manager connects the clusters accordingly.
 */
class FaceClusterFactory {
  public:
  FaceClusterFactory() = default;
  virtual ~FaceClusterFactory() = default;

  FaceClusterFactory(const FaceClusterFactory&) = delete;
  FaceClusterFactory(FaceClusterFactory&&) = delete;
  auto operator=(const FaceClusterFactory&) = delete;
  auto operator=(FaceClusterFactory&&) = delete;

  /// the faces of the layer `layerId`; a layer without faces gets no cluster
  [[nodiscard]] virtual std::size_t faceCount(std::size_t layerId) const = 0;

  /// before the clusters get made, with their time steps; on all ranks together
  virtual void prepare(const initializer::ClusterLayout& /*clusterLayout*/) {}

  /// the cluster of the faces of the layer `layerId`, which has faces; it takes the time steps of
  /// its time cluster in `clusterLayout`, and works on the executor of the cells next to it
  virtual std::unique_ptr<FaceCluster> makeCluster(std::size_t layerId,
                                                   const initializer::ClusterLayout& clusterLayout,
                                                   Executor executor,
                                                   const ClusterMonitoring& monitoring) = 0;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_FACECLUSTERFACTORY_H_
