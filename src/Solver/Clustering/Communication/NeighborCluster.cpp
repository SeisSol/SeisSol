// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "NeighborCluster.h"

#include <Initializer/BasicTypedefs.h>
#include <Parallel/Runtime/Stream.h>
namespace seissol::solver::clustering::communication {

HaloCommunication getHaloCommunication(const initializer::ClusterLayout& layout,
                                       const MeshStructure* structure) {
  HaloCommunication communication;
  communication.copy.resize(layout.globalClusterCount);
  communication.ghost.resize(layout.globalClusterCount);
  for (std::size_t i = 0; i < layout.localClusterIds.size(); ++i) {
    const auto& clusterStructure = structure[i];
    const std::size_t localClusterIndex = layout.localClusterIds[i];
    for (std::size_t j = 0; j < clusterStructure.numberOfRegions; ++j) {
      const std::size_t remoteClusterIndex = clusterStructure.neighboringClusters[j][1];
      // TODO: neighboring global-only clusters
      communication.copy.at(localClusterIndex)
          .emplace_back(RemoteCluster{clusterStructure.copyRegions[j],
                                      clusterStructure.copyRegionSizes[j],
                                      MPI_C_REAL,
                                      clusterStructure.neighboringClusters[j][0],
                                      DataTagOffset + clusterStructure.sendIdentifiers[j]});
      communication.ghost.at(remoteClusterIndex)
          .emplace_back(RemoteCluster{clusterStructure.ghostRegions[j],
                                      clusterStructure.ghostRegionSizes[j],
                                      MPI_C_REAL,
                                      clusterStructure.neighboringClusters[j][0],
                                      DataTagOffset + clusterStructure.receiveIdentifiers[j]});
    }
  }
  return communication;
}

NeighborCluster::NeighborCluster(const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                                 double priority)
    : myRuntime(parallel::runtime::StreamRuntime(cpuExecutor, priority)) {}

void NeighborCluster::startFrom(parallel::runtime::StreamRuntime& runtime) {
  void* event = runtime.recordEvent();
  this->myRuntime.waitEvent(event);
  start(this->myRuntime);
}

void NeighborCluster::stopTo(parallel::runtime::StreamRuntime& runtime) {
  stop(this->myRuntime);
  void* event = this->myRuntime.recordEvent();
  runtime.waitEvent(event);
}

} // namespace seissol::solver::clustering::communication
