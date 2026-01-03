// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifdef ACL_DEVICE

#include "Solver/TimeStepping/GhostTimeClusterWithCopy.h"

#include "Parallel/MPI.h"

#include <Device/device.h>

namespace seissol::time_stepping {
template <Mpi::DataTransferMode CommType>
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
  prefetchCopyRegionsStreams.resize(this->meshStructure_.copy.size());
  prefetchGhostRegionsStreams.resize(this->meshStructure_.ghost.size());
  receiveRegionsStates.resize(this->meshStructure_.ghost.size());
  duplicatedCopyRegions.resize(this->meshStructure_.copy.size());
  duplicatedGhostRegions.resize(this->meshStructure_.ghost.size());

  for (size_t region = 0; region < this->meshStructure_.copy.size(); ++region) {
    prefetchCopyRegionsStreams[region] = device.api->createStream();
    const size_t copyRegionSize = this->meshStructure_.copy[region].size *
                                  sizeOfRealType(this->meshStructure_.copy[region].datatype);
    if constexpr (CommType == Mpi::DataTransferMode::CopyInCopyOutHost) {
      duplicatedCopyRegions[region] = device.api->allocPinnedMem(copyRegionSize);
    }

    if (persistent) {
      MPI_Send_init(this->meshStructure_.copy[region].data,
                    static_cast<int>(this->meshStructure_.copy[region].size),
                    Mpi::precisionToMpiType(this->meshStructure_.copy[region].datatype),
                    this->meshStructure_.copy[region].rank,
                    this->meshStructure_.copy[region].tag,
                    seissol::Mpi::mpi.comm(),
                    sendRequests_.data() + region);
    }
  }
  for (size_t region = 0; region < this->meshStructure_.ghost.size(); ++region) {
    prefetchGhostRegionsStreams[region] = device.api->createStream();
    receiveRegionsStates[region] = ReceiveState::RequiresMpiTesting;
    const size_t ghostRegionSize = this->meshStructure_.ghost[region].size *
                                   sizeOfRealType(this->meshStructure_.ghost[region].datatype);
    if constexpr (CommType == Mpi::DataTransferMode::CopyInCopyOutHost) {
      duplicatedGhostRegions[region] = device.api->allocPinnedMem(ghostRegionSize);
    }

    if (persistent) {
      MPI_Recv_init(this->meshStructure_.ghost[region].data,
                    static_cast<int>(this->meshStructure_.ghost[region].size),
                    Mpi::precisionToMpiType(this->meshStructure_.ghost[region].datatype),
                    this->meshStructure_.ghost[region].rank,
                    this->meshStructure_.ghost[region].tag,
                    seissol::Mpi::mpi.comm(),
                    recvRequests_.data() + region);
    }
  }
}

template <Mpi::DataTransferMode CommType>
GhostTimeClusterWithCopy<CommType>::~GhostTimeClusterWithCopy() {
  for (size_t region = 0; region < this->meshStructure_.copy.size(); ++region) {
    device.api->destroyGenericStream(prefetchCopyRegionsStreams[region]);
    if constexpr (CommType == Mpi::DataTransferMode::CopyInCopyOutHost) {
      device.api->freePinnedMem(duplicatedCopyRegions[region]);
    }
  }
  for (size_t region = 0; region < this->meshStructure_.ghost.size(); ++region) {
    device.api->destroyGenericStream(prefetchGhostRegionsStreams[region]);
    if constexpr (CommType == Mpi::DataTransferMode::CopyInCopyOutHost) {
      device.api->freePinnedMem(duplicatedGhostRegions[region]);
    }
  }
}

template <Mpi::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::finalize() {
  if (persistent) {
    for (size_t region = 0; region < sendRequests_.size(); ++region) {
      MPI_Request_free(sendRequests_.data() + region);
    }
    for (size_t region = 0; region < recvRequests_.size(); ++region) {
      MPI_Request_free(recvRequests_.data() + region);
    }
  }
}

template <Mpi::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct_.correctionTime > lastSendTime_);

  lastSendTime_ = ct_.correctionTime;
  auto prefetchedRegions = prefetchCopyLayer();

  while (!prefetchedRegions.empty()) {
    for (auto region = prefetchedRegions.begin(); region != prefetchedRegions.end();) {
      auto* stream = prefetchCopyRegionsStreams[*region];
      if (device.api->isStreamWorkDone(stream)) {
        if (persistent) {
          MPI_Start(sendRequests_.data() + (*region));
        } else {
          MPI_Isend(this->meshStructure_.copy[*region].data,
                    static_cast<int>(this->meshStructure_.copy[*region].size),
                    Mpi::precisionToMpiType(this->meshStructure_.copy[*region].datatype),
                    this->meshStructure_.copy[*region].rank,
                    this->meshStructure_.copy[*region].tag,
                    seissol::Mpi::mpi.comm(),
                    sendRequests_.data() + (*region));
        }
        sendQueue_.push_back(*region);
        region = prefetchedRegions.erase(region);
      } else {
        ++region;
      }
    }
  }
}

template <Mpi::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct_.predictionTime >= lastSendTime_);
  if (persistent) {
    MPI_Startall(recvRequests_.size(), recvRequests_.data());
  }
  for (std::size_t region = 0; region < recvRequests_.size(); ++region) {
    if (!persistent) {
      MPI_Irecv(this->meshStructure_.ghost[region].data,
                static_cast<int>(this->meshStructure_.ghost[region].size),
                Mpi::precisionToMpiType(this->meshStructure_.ghost[region].datatype),
                this->meshStructure_.ghost[region].rank,
                this->meshStructure_.ghost[region].tag,
                seissol::Mpi::mpi.comm(),
                recvRequests_.data() + region);
    }
    receiveRegionsStates[region] = ReceiveState::RequiresMpiTesting;
    receiveQueue_.push_back(region);
  }
}

template <Mpi::DataTransferMode CommType>
bool GhostTimeClusterWithCopy<CommType>::testReceiveQueue() {
  for (auto region = receiveQueue_.begin(); region != receiveQueue_.end();) {
    const auto state = receiveRegionsStates[*region];

    switch (state) {
    case ReceiveState::RequiresMpiTesting: {
      int testSuccess = 0;
      MPI_Request* request = &(recvRequests_.data())[*region];
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
        region = receiveQueue_.erase(region);
      } else {
        ++region;
      }
      break;
    }
    case ReceiveState::Ready: {
      region = receiveQueue_.erase(region);
    }
    }
  }
  return receiveQueue_.empty();
}

template <Mpi::DataTransferMode CommType>
bool GhostTimeClusterWithCopy<CommType>::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testReceiveQueue();
}

template <Mpi::DataTransferMode CommType>
std::list<int> GhostTimeClusterWithCopy<CommType>::prefetchCopyLayer() {
  std::list<int> prefetchedRegions{};
  for (std::size_t region = 0; region < prefetchCopyRegionsStreams.size(); ++region) {
    auto* stream = prefetchCopyRegionsStreams[region];
    const auto messageSize = this->meshStructure_.copy[region].size;

    if constexpr (CommType == Mpi::DataTransferMode::CopyInCopyOutHost) {
      device.api->copyFromAsync(duplicatedCopyRegions[region],
                                this->meshStructure_.copy[region].data,
                                messageSize *
                                    sizeOfRealType(this->meshStructure_.copy[region].datatype),
                                stream);
    }
    prefetchedRegions.push_back(region);
  }
  return prefetchedRegions;
}

template <Mpi::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::prefetchGhostRegion(std::size_t region) {
  auto* stream = prefetchGhostRegionsStreams[region];
  const auto messageSize = this->meshStructure_.ghost[region].size;
  if constexpr (CommType == Mpi::DataTransferMode::CopyInCopyOutHost) {
    device.api->copyToAsync(this->meshStructure_.ghost[region].data,
                            duplicatedGhostRegions[region],
                            messageSize *
                                sizeOfRealType(this->meshStructure_.ghost[region].datatype),
                            stream);
  }
}

template class GhostTimeClusterWithCopy<Mpi::DataTransferMode::CopyInCopyOutHost>;
} // namespace seissol::time_stepping
#endif // ACL_DEVICE
