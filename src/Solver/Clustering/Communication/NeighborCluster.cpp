// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "NeighborCluster.h"

#include <Initializer/BasicTypedefs.h>
#include <Parallel/Runtime/Stream.h>
namespace seissol::solver::clustering::communication {

HaloCommunication getHaloCommunication(std::size_t clusterCount, const MeshStructure* structure) {
  HaloCommunication communication;
  communication.copy.resize(clusterCount);
  communication.ghost.resize(clusterCount);
  for (std::size_t i = 0; i < clusterCount; ++i) {
    const auto& clusterStructure = structure[i];
    for (std::size_t j = 0; j < clusterStructure.numberOfRegions; ++j) {
      const std::size_t clusterIndex = clusterStructure.neighboringClusters[j][0];
      communication.copy.at(clusterIndex)
          .emplace_back(RemoteCluster{clusterStructure.copyRegions[j],
                                      clusterStructure.copyRegionSizes[j],
                                      MPI_C_REAL,
                                      clusterStructure.neighboringClusters[j][0],
                                      DataTagOffset + clusterStructure.sendIdentifiers[j]});
      communication.ghost.at(clusterIndex)
          .emplace_back(RemoteCluster{clusterStructure.ghostRegions[j],
                                      clusterStructure.ghostRegionSizes[j],
                                      MPI_C_REAL,
                                      clusterStructure.neighboringClusters[j][0],
                                      DataTagOffset + clusterStructure.receiveIdentifiers[j]});
    }
  }
  return communication;
}

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
