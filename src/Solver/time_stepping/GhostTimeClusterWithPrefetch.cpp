#include "Initializer/typedefs.hpp"
#ifdef REQUIRED_COMM_LAYERS_PREFETCH

#include <Parallel/MPI.h>
#include "Solver/time_stepping/GhostTimeClusterWithPrefetch.h"
#include "device.h"

namespace seissol::time_stepping {
void GhostTimeClusterWithPrefetch::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.correctionTime > lastSendTime);

  lastSendTime = ct.correctionTime;
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      syncCommRegion(region, Direction::Send);
      MPI_Isend(meshStructure->duplicatedCopyRegions[region],
                static_cast<int>(meshStructure->copyRegionSizes[region]),
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->sendIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->sendRequests + region);
      sendQueue.push_back(region);
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
      receiveQueue.push_back(region);
    }
  }
}

bool GhostTimeClusterWithPrefetch::testReceiveQueue() {
  for (auto region = receiveQueue.begin(); region != receiveQueue.end();) {
    MPI_Request* request = &(meshStructure->receiveRequests)[*region];
    int testSuccess = 0;
    MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess) {
      syncCommRegion(*region, Direction::Recv);
      region = receiveQueue.erase(region);
    } else {
      ++region;
    }
  }
  return receiveQueue.empty();
}

bool GhostTimeClusterWithPrefetch::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testReceiveQueue();
}

void GhostTimeClusterWithPrefetch::syncCommRegion(int region, Direction direction) {
  device::DeviceInstance& device = device::DeviceInstance::getInstance();

  switch (direction) {
  case Direction::Send: {
    const auto messageSize = static_cast<int>(meshStructure->copyRegionSizes[region]);
#if defined(PREFETCH_COMM_LAYERS_TO_HOST)
    device.api->copyFrom(meshStructure->duplicatedCopyRegions[region],
                         meshStructure->copyRegions[region],
                         messageSize * sizeof(real));
#elif defined(PREFETCH_COMM_LAYERS_TO_DEVICE)
    device.api->copyBetween(meshStructure->duplicatedCopyRegions[region],
                            meshStructure->copyRegions[region],
                            messageSize * sizeof(real));
#else
#error "comm layers prefetch type is not selected"
#endif
    break;
  }
  case Direction::Recv: {
    const auto messageSize = static_cast<int>(meshStructure->ghostRegionSizes[region]);
#if defined(PREFETCH_COMM_LAYERS_TO_HOST)
    device.api->copyTo(meshStructure->ghostRegions[region],
                       meshStructure->duplicatedGhostRegions[region],
                       messageSize * sizeof(real));
#elif defined(PREFETCH_COMM_LAYERS_TO_DEVICE)
    device.api->copyBetween(meshStructure->ghostRegions[region],
                            meshStructure->duplicatedGhostRegions[region],
                            messageSize * sizeof(real));
#else
#error "comm layers prefetch type is not selected"
#endif
    break;
  }
  }
}
} // namespace seissol::time_stepping
#endif // REQUIRED_COMM_LAYERS_PREFETCH
