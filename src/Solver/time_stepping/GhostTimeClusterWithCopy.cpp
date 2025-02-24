// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Typedefs.h"
#ifdef ACL_DEVICE

#include "Parallel/MPI.h"
#include "Solver/time_stepping/GhostTimeClusterWithCopy.h"
#include "device.h"

namespace seissol::time_stepping {
template <MPI::DataTransferMode CommType>
GhostTimeClusterWithCopy<CommType>::GhostTimeClusterWithCopy(
    double maxTimeStepSize,
    int timeStepRate,
    int globalTimeClusterId,
    int otherGlobalTimeClusterId,
    const MeshStructure* meshStructure,
    bool persistent)
    : AbstractGhostTimeCluster(maxTimeStepSize,
                               timeStepRate,
                               globalTimeClusterId,
                               otherGlobalTimeClusterId,
                               meshStructure), persistent(persistent) {
  numberOfRegions = meshStructure->numberOfRegions;
  prefetchCopyRegionsStreams.resize(numberOfRegions);
  prefetchGhostRegionsStreams.resize(numberOfRegions);
  receiveRegionsStates.resize(numberOfRegions);

  for (size_t region = 0; region < numberOfRegions; ++region) {
    prefetchCopyRegionsStreams[region] = device.api->createStream();
    prefetchGhostRegionsStreams[region] = device.api->createStream();
    receiveRegionsStates[region] = ReceiveState::RequiresMpiTesting;
  }

  duplicatedCopyRegions.resize(numberOfRegions);
  duplicatedGhostRegions.resize(numberOfRegions);
  for (size_t region = 0; region < numberOfRegions; ++region) {
    const size_t copyRegionSize = meshStructure->copyRegionSizes[region] * sizeof(real);
    const size_t ghostRegionSize = meshStructure->ghostRegionSizes[region] * sizeof(real);
    if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
      duplicatedCopyRegions[region] =
          static_cast<real*>(device.api->allocPinnedMem(copyRegionSize));
      duplicatedGhostRegions[region] =
          static_cast<real*>(device.api->allocPinnedMem(ghostRegionSize));
    }

    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      if (persistent) {
        MPI_Send_init(duplicatedCopyRegions[region],
                  static_cast<int>(meshStructure->copyRegionSizes[region]),
                  MPI_C_REAL,
                  meshStructure->neighboringClusters[region][0],
                  DataTagOffset + meshStructure->sendIdentifiers[region],
                  seissol::MPI::mpi.comm(),
                  meshStructure->sendRequests + region);
        MPI_Recv_init(duplicatedGhostRegions[region],
                  static_cast<int>(meshStructure->ghostRegionSizes[region]),
                  MPI_C_REAL,
                  meshStructure->neighboringClusters[region][0],
                  DataTagOffset + meshStructure->receiveIdentifiers[region],
                  seissol::MPI::mpi.comm(),
                  meshStructure->receiveRequests + region);
      }
    }
  }
}

template <MPI::DataTransferMode CommType>
GhostTimeClusterWithCopy<CommType>::~GhostTimeClusterWithCopy() {
  for (size_t region = 0; region < numberOfRegions; ++region) {
    device.api->destroyGenericStream(prefetchCopyRegionsStreams[region]);
    device.api->destroyGenericStream(prefetchGhostRegionsStreams[region]);
    if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
      device.api->freePinnedMem(duplicatedCopyRegions[region]);
      device.api->freePinnedMem(duplicatedGhostRegions[region]);
    }
  }
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::finalize() {
  if (persistent) {
    for (size_t region = 0; region < numberOfRegions; ++region) {
      if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
        MPI_Request_free(meshStructure->sendRequests + region);
        MPI_Request_free(meshStructure->receiveRequests + region);
      }
    }
  }
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.correctionTime > lastSendTime);

  lastSendTime = ct.correctionTime;
  auto prefetchedRegions = prefetchCopyLayer();

  while (!prefetchedRegions.empty()) {
    for (auto region = prefetchedRegions.begin(); region != prefetchedRegions.end();) {
      auto* stream = prefetchCopyRegionsStreams[*region];
      if (device.api->isStreamWorkDone(stream)) {
        if (persistent) {
          MPI_Start(meshStructure->sendRequests + (*region));
        }
        else {
          MPI_Isend(duplicatedCopyRegions[*region],
                    static_cast<int>(meshStructure->copyRegionSizes[*region]),
                    MPI_C_REAL,
                    meshStructure->neighboringClusters[*region][0],
                    DataTagOffset + meshStructure->sendIdentifiers[*region],
                    seissol::MPI::mpi.comm(),
                    meshStructure->sendRequests + (*region));
        }
        sendQueue.push_back(*region);
        region = prefetchedRegions.erase(region);
      } else {
        ++region;
      }
    }
  }
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.predictionTime >= lastSendTime);
  for (unsigned int region = 0; region < numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      if (persistent) {
        MPI_Start(meshStructure->receiveRequests + region);
      }
      else {
        MPI_Irecv(duplicatedGhostRegions[region],
                  static_cast<int>(meshStructure->ghostRegionSizes[region]),
                  MPI_C_REAL,
                  meshStructure->neighboringClusters[region][0],
                  DataTagOffset + meshStructure->receiveIdentifiers[region],
                  seissol::MPI::mpi.comm(),
                  meshStructure->receiveRequests + region);
      }
      receiveRegionsStates[region] = ReceiveState::RequiresMpiTesting;
      receiveQueue.push_back(region);
    }
  }
}

template <MPI::DataTransferMode CommType>
bool GhostTimeClusterWithCopy<CommType>::testReceiveQueue() {
  for (auto region = receiveQueue.begin(); region != receiveQueue.end();) {
    const auto state = receiveRegionsStates[*region];

    switch (state) {
    case ReceiveState::RequiresMpiTesting: {
      int testSuccess = 0;
      MPI_Request* request = &(meshStructure->receiveRequests)[*region];
      MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
      if (testSuccess) {
        prefetchGhostRegion(*region);
        receiveRegionsStates[*region] = ReceiveState::RequiresPrefetchTesting;
      }
      ++region;
      break;
    }
    case ReceiveState::RequiresPrefetchTesting: {
      auto* stream = prefetchGhostRegionsStreams[*region];
      if (device.api->isStreamWorkDone(stream)) {
        receiveRegionsStates[*region] = ReceiveState::Ready;
        // Note: fall-through to the `Ready` state is intentional
      } else {
        ++region;
        break;
      }
    }
    case ReceiveState::Ready: {
      region = receiveQueue.erase(region);
    }
    }
  }
  return receiveQueue.empty();
}

template <MPI::DataTransferMode CommType>
bool GhostTimeClusterWithCopy<CommType>::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testReceiveQueue();
}

template <MPI::DataTransferMode CommType>
std::list<int> GhostTimeClusterWithCopy<CommType>::prefetchCopyLayer() {
  std::list<int> prefetchedRegions{};
  for (unsigned int region = 0; region < numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {

      auto* stream = prefetchCopyRegionsStreams[region];
      const auto messageSize = static_cast<int>(meshStructure->copyRegionSizes[region]);

      if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
        device.api->copyFromAsync(duplicatedCopyRegions[region],
                                  meshStructure->copyRegions[region],
                                  messageSize * sizeof(real),
                                  stream);
      }
      prefetchedRegions.push_back(region);
    }
  }
  return prefetchedRegions;
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::prefetchGhostRegion(int region) {
  auto* stream = prefetchGhostRegionsStreams[region];
  const auto messageSize = static_cast<int>(meshStructure->ghostRegionSizes[region]);
  if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
    device.api->copyToAsync(meshStructure->ghostRegions[region],
                            duplicatedGhostRegions[region],
                            messageSize * sizeof(real),
                            stream);
  }
}

template class GhostTimeClusterWithCopy<MPI::DataTransferMode::CopyInCopyOutHost>;
} // namespace seissol::time_stepping
#endif // ACL_DEVICE

