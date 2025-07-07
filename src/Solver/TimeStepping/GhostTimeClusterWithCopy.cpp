// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifdef ACL_DEVICE

#include "Solver/TimeStepping/GhostTimeClusterWithCopy.h"
#include "Parallel/MPI.h"
#include "device.h"

namespace seissol::time_stepping {
template <MPI::DataTransferMode CommType>
GhostTimeClusterWithCopy<CommType>::GhostTimeClusterWithCopy(
    double maxTimeStepSize,
    int timeStepRate,
    int globalTimeClusterId,
    int otherGlobalTimeClusterId,
    const seissol::solver::HaloCommunication& meshStructure,
    bool persistent)
    : AbstractGhostTimeCluster(maxTimeStepSize,
                               timeStepRate,
                               globalTimeClusterId,
                               otherGlobalTimeClusterId,
                               meshStructure),
      persistent(persistent) {
  numberOfRegions = this->meshStructure.size();
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
    const size_t copyRegionSize = this->meshStructure[region].copy.size *
                                  sizeOfRealType(this->meshStructure[region].copy.datatype);
    const size_t ghostRegionSize = this->meshStructure[region].ghost.size *
                                   sizeOfRealType(this->meshStructure[region].ghost.datatype);
    if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
      duplicatedCopyRegions[region] = device.api->allocPinnedMem(copyRegionSize);
      duplicatedGhostRegions[region] = device.api->allocPinnedMem(ghostRegionSize);
    }

    if (persistent) {
      MPI_Send_init(this->meshStructure[region].copy.data,
                    static_cast<int>(this->meshStructure[region].copy.size),
                    MPI::precisionToMpiType(this->meshStructure[region].copy.datatype),
                    this->meshStructure[region].copy.rank,
                    this->meshStructure[region].copy.tag,
                    seissol::MPI::mpi.comm(),
                    sendRequests.data() + region);
      MPI_Recv_init(this->meshStructure[region].ghost.data,
                    static_cast<int>(this->meshStructure[region].ghost.size),
                    MPI::precisionToMpiType(this->meshStructure[region].ghost.datatype),
                    this->meshStructure[region].ghost.rank,
                    this->meshStructure[region].ghost.tag,
                    seissol::MPI::mpi.comm(),
                    recvRequests.data() + region);
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
      MPI_Request_free(sendRequests.data() + region);
      MPI_Request_free(recvRequests.data() + region);
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
          MPI_Start(sendRequests.data() + (*region));
        } else {
          MPI_Isend(this->meshStructure[*region].copy.data,
                    static_cast<int>(this->meshStructure[*region].copy.size),
                    MPI::precisionToMpiType(this->meshStructure[*region].copy.datatype),
                    this->meshStructure[*region].copy.rank,
                    this->meshStructure[*region].copy.tag,
                    seissol::MPI::mpi.comm(),
                    sendRequests.data() + (*region));
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
  if (persistent) {
    MPI_Startall(recvRequests.size(), recvRequests.data());
  }
  for (std::size_t region = 0; region < numberOfRegions; ++region) {
    if (!persistent) {
      MPI_Irecv(this->meshStructure[region].ghost.data,
                static_cast<int>(this->meshStructure[region].ghost.size),
                MPI::precisionToMpiType(this->meshStructure[region].ghost.datatype),
                this->meshStructure[region].ghost.rank,
                this->meshStructure[region].ghost.tag,
                seissol::MPI::mpi.comm(),
                recvRequests.data() + region);
    }
    receiveRegionsStates[region] = ReceiveState::RequiresMpiTesting;
    receiveQueue.push_back(region);
  }
}

template <MPI::DataTransferMode CommType>
bool GhostTimeClusterWithCopy<CommType>::testReceiveQueue() {
  for (auto region = receiveQueue.begin(); region != receiveQueue.end();) {
    const auto state = receiveRegionsStates[*region];

    switch (state) {
    case ReceiveState::RequiresMpiTesting: {
      int testSuccess = 0;
      MPI_Request* request = &(recvRequests.data())[*region];
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
  for (std::size_t region = 0; region < numberOfRegions; ++region) {
    auto* stream = prefetchCopyRegionsStreams[region];
    const auto messageSize = this->meshStructure[region].copy.size;

    if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
      device.api->copyFromAsync(duplicatedCopyRegions[region],
                                this->meshStructure[region].copy.data,
                                messageSize *
                                    sizeOfRealType(this->meshStructure[region].copy.datatype),
                                stream);
    }
    prefetchedRegions.push_back(region);
  }
  return prefetchedRegions;
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::prefetchGhostRegion(std::size_t region) {
  auto* stream = prefetchGhostRegionsStreams[region];
  const auto messageSize = this->meshStructure[region].ghost.size;
  if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
    device.api->copyToAsync(this->meshStructure[region].ghost.data,
                            duplicatedGhostRegions[region],
                            messageSize *
                                sizeOfRealType(this->meshStructure[region].ghost.datatype),
                            stream);
  }
}

template class GhostTimeClusterWithCopy<MPI::DataTransferMode::CopyInCopyOutHost>;
} // namespace seissol::time_stepping
#endif // ACL_DEVICE
