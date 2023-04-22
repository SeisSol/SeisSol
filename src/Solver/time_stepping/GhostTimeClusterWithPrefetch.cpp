#include "Initializer/typedefs.hpp"
#ifdef REQUIRED_COMM_LAYERS_PREFETCH

#include <Parallel/MPI.h>
#include "Solver/time_stepping/GhostTimeClusterWithPrefetch.h"
#include "device.h"


namespace seissol::time_stepping {
GhostTimeClusterWithPrefetch::GhostTimeClusterWithPrefetch(double maxTimeStepSize,
                                                           int timeStepRate,
                                                           int globalTimeClusterId,
                                                           int otherGlobalTimeClusterId,
                                                           const MeshStructure* meshStructure)
    : GenericGhostTimeCluster(maxTimeStepSize,
                              timeStepRate,
                              globalTimeClusterId,
                              otherGlobalTimeClusterId,
                              meshStructure) {
  prefetchCopyRegionsStreams.resize(meshStructure->numberOfRegions);
  prefetchGhostRegionsStreams.resize(meshStructure->numberOfRegions);
  receiveRegionsStates.resize(meshStructure->numberOfRegions);

  for (size_t region = 0; region < meshStructure->numberOfRegions; ++region) {
    prefetchCopyRegionsStreams[region] = device.api->createGenericStream();
    prefetchGhostRegionsStreams[region] = device.api->createGenericStream();
    receiveRegionsStates[region] = RequiresMpiTesting;
  }
}

GhostTimeClusterWithPrefetch::~GhostTimeClusterWithPrefetch() {
  for (size_t region = 0; region < meshStructure->numberOfRegions; ++region) {
    device.api->destroyGenericStream(prefetchCopyRegionsStreams[region]);
    device.api->destroyGenericStream(prefetchGhostRegionsStreams[region]);
  }
}

void GhostTimeClusterWithPrefetch::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.correctionTime > lastSendTime);

  lastSendTime = ct.correctionTime;
  auto prefetchedRegions = prefetchCopyLayer();

  while (not prefetchedRegions.empty()) {
    for (auto region = prefetchedRegions.begin(); region != prefetchedRegions.end();) {
      auto* stream = prefetchCopyRegionsStreams[*region];
      if (device.api->isStreamWorkDone(stream)) {
        MPI_Isend(meshStructure->duplicatedCopyRegions[*region],
                  static_cast<int>(meshStructure->copyRegionSizes[*region]),
                  MPI_C_REAL,
                  meshStructure->neighboringClusters[*region][0],
                  timeData + meshStructure->sendIdentifiers[*region],
                  seissol::MPI::mpi.comm(),
                  meshStructure->sendRequests + (*region));
        sendQueue.push_back(*region);
        region = prefetchedRegions.erase(region);
      } else {
        ++region;
      }
    }
  }
}

void GhostTimeClusterWithPrefetch::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.predictionTime > lastSendTime);
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      MPI_Irecv(meshStructure->duplicatedGhostRegions[region],
                static_cast<int>(meshStructure->ghostRegionSizes[region]),
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->receiveIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->receiveRequests + region);
      receiveRegionsStates[region] = RequiresMpiTesting;
      receiveQueue.push_back(region);
    }
  }
}

bool GhostTimeClusterWithPrefetch::testReceiveQueue() {
  for (auto region = receiveQueue.begin(); region != receiveQueue.end();) {
    const auto state = receiveRegionsStates[*region];

    switch(state) {
    case ReceiveState::RequiresMpiTesting:  {
      int testSuccess = 0;
      MPI_Request* request = &(meshStructure->receiveRequests)[*region];
      MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
      if (testSuccess) {
        prefetchGhostRegion(*region);
        receiveRegionsStates[*region] = RequiresPrefetchTesting;
      }
      ++region;
      break;
    }
    case ReceiveState::RequiresPrefetchTesting: {
      auto* stream = prefetchGhostRegionsStreams[*region];
      if (device.api->isStreamWorkDone(stream)) {
        receiveRegionsStates[*region] = Ready;
        // Note: fall to `Ready` is intentional
      }
      else {
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

bool GhostTimeClusterWithPrefetch::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testReceiveQueue();
}

std::list<int> GhostTimeClusterWithPrefetch::prefetchCopyLayer() {
  std::list<int> prefetchedRegions{};
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {

      auto* stream = prefetchCopyRegionsStreams[region];
      const auto messageSize = static_cast<int>(meshStructure->copyRegionSizes[region]);

#if defined(PREFETCH_COMM_LAYERS_TO_HOST)
      device.api->copyFromAsync(meshStructure->duplicatedCopyRegions[region],
                                meshStructure->copyRegions[region],
                                messageSize * sizeof(real),
                                stream);
#elif defined(PREFETCH_COMM_LAYERS_TO_DEVICE)
      device.api->copyBetweenAsync(meshStructure->duplicatedCopyRegions[region],
                                   meshStructure->copyRegions[region],
                                   messageSize * sizeof(real),
                                   stream);
#else
#error "comm layers prefetch type is not selected"
#endif

      prefetchedRegions.push_back(region);
    }
  }
  return prefetchedRegions;
}

void GhostTimeClusterWithPrefetch::prefetchGhostRegion(int region) {
  auto* stream = prefetchGhostRegionsStreams[region];
  const auto messageSize = static_cast<int>(meshStructure->ghostRegionSizes[region]);
#if defined(PREFETCH_COMM_LAYERS_TO_HOST)
  device.api->copyToAsync(meshStructure->ghostRegions[region],
                          meshStructure->duplicatedGhostRegions[region],
                          messageSize * sizeof(real),
                          stream);
#elif defined(PREFETCH_COMM_LAYERS_TO_DEVICE)
  device.api->copyBetweenAsync(meshStructure->ghostRegions[region],
                               meshStructure->duplicatedGhostRegions[region],
                               messageSize * sizeof(real),
                               stream);
#else
#error "comm layers prefetch type is not selected"
#endif
}
} // namespace seissol::time_stepping
#endif // REQUIRED_COMM_LAYERS_PREFETCH
