#include "Initializer/typedefs.hpp"
#include <Initializer/BasicTypedefs.hpp>
#include <Initializer/preProcessorMacros.hpp>
#include <Kernels/precision.hpp>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/AbstractGhostTimeCluster.h>
#include <cassert>
#include <cstddef>
#include <list>
#include <mpi.h>
#ifdef ACL_DEVICE

#include "Parallel/MPI.h"
#include "Solver/Clustering/Communication/GhostTimeClusterWithCopy.h"
#include "device.h"

namespace seissol::time_stepping {
template <MPI::DataTransferMode CommType>
GhostTimeClusterWithCopy<CommType>::GhostTimeClusterWithCopy(double maxTimeStepSize,
                                                             int timeStepRate,
                                                             int globalTimeClusterId,
                                                             int otherGlobalTimeClusterId,
                                                             const MeshStructure* meshStructure,
                                                             bool persistent)
    : AbstractGhostTimeCluster(maxTimeStepSize,
                               timeStepRate,
                               globalTimeClusterId,
                               otherGlobalTimeClusterId,
                               meshStructure),
      persistent(persistent) {

  const std::size_t neighbors = ghostClusters.size();

  prefetchCopyRegionsStreams.resize(neighbors);
  prefetchGhostRegionsStreams.resize(neighbors);
  receiveRegionsStates.resize(neighbors);

  for (std::size_t i = 0; i < neighbors; ++i) {
    prefetchCopyRegionsStreams[i] = device.api->createGenericStream();
    prefetchGhostRegionsStreams[i] = device.api->createGenericStream();
    receiveRegionsStates[i] = ReceiveState::RequiresMpiTesting;
  }

  duplicatedCopyRegions.resize(neighbors);
  duplicatedGhostRegions.resize(neighbors);
  for (std::size_t i = 0; i < neighbors; ++i) {
    const size_t copyRegionSize = copyClusters[i].size * sizeof(real);
    const size_t ghostRegionSize = ghostClusters[i].size * sizeof(real);
    if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
      duplicatedCopyRegions[i] = static_cast<real*>(device.api->allocPinnedMem(copyRegionSize));
      duplicatedGhostRegions[i] = static_cast<real*>(device.api->allocPinnedMem(ghostRegionSize));
    } else {
      static_assert(false, "Unhandled MPI copy case.");
    }

    if (persistent) {
      {
        const auto& cluster = copyClusters.at(i);
        MPI_Send_init(duplicatedCopyRegions[i],
                      static_cast<int>(cluster.size),
                      cluster.datatype,
                      cluster.rank,
                      cluster.tag,
                      seissol::MPI::mpi.comm(),
                      &sendRequests[i]);
      }
      {
        const auto& cluster = ghostClusters.at(i);
        MPI_Recv_init(duplicatedGhostRegions[i],
                      static_cast<int>(cluster.size),
                      cluster.datatype,
                      cluster.rank,
                      cluster.tag,
                      seissol::MPI::mpi.comm(),
                      &recvRequests[i]);
      }
    }
  }
}

template <MPI::DataTransferMode CommType>
GhostTimeClusterWithCopy<CommType>::~GhostTimeClusterWithCopy() {
  const std::size_t neighbors = ghostClusters.size();
  for (std::size_t i = 0; i < neighbors; ++i) {
    device.api->destroyGenericStream(prefetchCopyRegionsStreams[i]);
    device.api->destroyGenericStream(prefetchGhostRegionsStreams[i]);
    if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
      device.api->freePinnedMem(duplicatedCopyRegions[i]);
      device.api->freePinnedMem(duplicatedGhostRegions[i]);
    }
  }
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::finalize() {
  if (persistent) {
    for (auto& request : sendRequests) {
      MPI_Request_free(&request);
    }
    for (auto& request : recvRequests) {
      MPI_Request_free(&request);
    }
  }
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.time.at(ComputeStep::Correct) > lastSendTime);

  lastSendTime = ct.time.at(ComputeStep::Correct);
  auto prefetchedRegions = prefetchCopyLayer();

  while (!prefetchedRegions.empty()) {
    for (auto it = prefetchedRegions.begin(); it != prefetchedRegions.end();) {
      auto* stream = prefetchCopyRegionsStreams[*it];
      if (device.api->isStreamWorkDone(stream)) {
        if (persistent) {
          MPI_Start(&sendRequests[*it]);
        } else {
          const auto& cluster = copyClusters[*it];
          MPI_Isend(duplicatedCopyRegions[*it],
                    static_cast<int>(cluster.size),
                    cluster.datatype,
                    cluster.rank,
                    cluster.tag,
                    seissol::MPI::mpi.comm(),
                    &sendRequests[*it]);
        }
        sendQueue.push_back(*it);
        it = prefetchedRegions.erase(it);
      } else {
        ++it;
      }
    }
  }
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.time.at(ComputeStep::Predict) >= lastSendTime);
  if (persistent) {
    MPI_Startall(recvRequests.size(), recvRequests.data());
    for (std::size_t i = 0; i < ghostClusters.size(); ++i) {
      recvQueue.push_back(i);
    }
  }
  for (std::size_t i = 0; i < ghostClusters.size(); ++i) {
    if (!persistent) {
      const auto& cluster = ghostClusters[i];
      MPI_Irecv(duplicatedGhostRegions[i],
                static_cast<int>(cluster.size),
                cluster.datatype,
                cluster.rank,
                cluster.tag,
                seissol::MPI::mpi.comm(),
                &recvRequests[i]);
    }
    receiveRegionsStates[i] = ReceiveState::RequiresMpiTesting;
    recvQueue.push_back(i);
  }
}

template <MPI::DataTransferMode CommType>
bool GhostTimeClusterWithCopy<CommType>::testReceiveQueue() {
  for (auto it = recvQueue.begin(); it != recvQueue.end();) {
    const auto state = receiveRegionsStates[*it];

    switch (state) {
    case ReceiveState::RequiresMpiTesting: {
      int testSuccess = 0;
      MPI_Request* request = &recvRequests[*it];
      MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
      if (testSuccess) {
        prefetchGhostRegion(*it);
        receiveRegionsStates[*it] = ReceiveState::RequiresPrefetchTesting;
      }
      ++it;
      break;
    }
    case ReceiveState::RequiresPrefetchTesting: {
      auto* stream = prefetchGhostRegionsStreams[*it];
      if (device.api->isStreamWorkDone(stream)) {
        receiveRegionsStates[*it] = ReceiveState::Ready;
        // Note: fall-through to the `Ready` state is intentional
      } else {
        ++it;
        break;
      }
    }
    case ReceiveState::Ready: {
      it = recvQueue.erase(it);
    }
    }
  }
  return recvQueue.empty();
}

template <MPI::DataTransferMode CommType>
bool GhostTimeClusterWithCopy<CommType>::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testReceiveQueue();
}

template <MPI::DataTransferMode CommType>
std::list<std::size_t> GhostTimeClusterWithCopy<CommType>::prefetchCopyLayer() {
  std::list<std::size_t> prefetchedRegions{};
  for (std::size_t i = 0; i < copyClusters.size(); ++i) {
    auto* stream = prefetchCopyRegionsStreams[i];
    const auto messageSize = static_cast<int>(copyClusters[i].size);

    if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
      device.api->copyFromAsync(
          duplicatedCopyRegions[i], copyClusters[i].data, messageSize * sizeof(real), stream);
    }
    prefetchedRegions.push_back(i);
  }
  return prefetchedRegions;
}

template <MPI::DataTransferMode CommType>
void GhostTimeClusterWithCopy<CommType>::prefetchGhostRegion(std::size_t i) {
  auto* stream = prefetchGhostRegionsStreams[i];
  const auto messageSize = static_cast<int>(ghostClusters[i].size);
  if constexpr (CommType == MPI::DataTransferMode::CopyInCopyOutHost) {
    device.api->copyToAsync(
        copyClusters[i].data, duplicatedGhostRegions[i], messageSize * sizeof(real), stream);
  }
}

template class GhostTimeClusterWithCopy<MPI::DataTransferMode::CopyInCopyOutHost>;
} // namespace seissol::time_stepping
#endif // ACL_DEVICE
