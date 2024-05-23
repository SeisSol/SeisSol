#include "AbstractGhostTimeCluster.h"
#include "Parallel/MPI.h"
#include "Solver/time_stepping/DirectGhostTimeCluster.h"

#ifdef USE_CCL
#include <rccl.h>

#if REAL_SIZE == 8
using CCLReal = ncclFloat64;
#elif REAL_SIZE == 4
using CCLReal = ncclFloat32;
#endif
#endif

namespace seissol::time_stepping {
void DirectGhostTimeCluster::sendCopyLayer() {
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  assert(ct.correctionTime > lastSendTime);
  lastSendTime = ct.correctionTime;
#ifdef USE_CCL
  ncclGroupStart();
#endif
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
#ifdef USE_CCL
        ncclSend(meshStructure->copyRegions[region],
          static_cast<size_t>(meshStructure->copyRegionSizes[region]),
          CCLReal,
          meshStructure->neighboringClusters[region][0],
          timeData + meshStructure->sendIdentifiers[region],
          static_cast<ncclComm_t>(comm),
          stream);
#else
      if (persistent) {
        MPI_Start(meshStructure->sendRequests + region);
      }
      else {
        MPI_Isend(meshStructure->copyRegions[region],
                    static_cast<int>(meshStructure->copyRegionSizes[region]),
                    MPI_C_REAL,
                    meshStructure->neighboringClusters[region][0],
                    timeData + meshStructure->sendIdentifiers[region],
                    seissol::MPI::mpi.comm(),
                    meshStructure->sendRequests + region
                  );
      }
      sendQueue.push_back(region);
#endif
    }
  }
#ifdef USE_CCL
  ncclGroupEnd();
#endif
}

void DirectGhostTimeCluster::receiveGhostLayer() {
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  assert(ct.predictionTime >= lastSendTime);
#ifdef USE_CCL
  ncclGroupStart();
#endif
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId) ) {
#ifdef USE_CCL
        ncclRecv(meshStructure->copyRegions[region],
          static_cast<size_t>(meshStructure->copyRegionSizes[region]),
          CCLReal,
          meshStructure->neighboringClusters[region][0],
          timeData + meshStructure->sendIdentifiers[region],
          static_cast<ncclComm_t>(comm),
          stream);
#else
      if (persistent) {
        MPI_Start(meshStructure->receiveRequests + region);
      }
      else {
        MPI_Irecv(meshStructure->ghostRegions[region],
                  static_cast<int>(meshStructure->ghostRegionSizes[region]),
                  MPI_C_REAL,
                  meshStructure->neighboringClusters[region][0],
                  timeData + meshStructure->receiveIdentifiers[region],
                  seissol::MPI::mpi.comm(),
                  meshStructure->receiveRequests + region);
      }
      receiveQueue.push_back(region);
#endif
    }
  }
#ifdef USE_CCL
  ncclGroupEnd();
#endif
}

bool DirectGhostTimeCluster::testForGhostLayerReceives() {
  SCOREP_USER_REGION( "testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION )
  return testQueue(meshStructure->receiveRequests, receiveQueue);
}

DirectGhostTimeCluster::DirectGhostTimeCluster(double maxTimeStepSize,
                                               int timeStepRate,
                                               int globalTimeClusterId,
                                               int otherGlobalTimeClusterId,
                                               const MeshStructure *meshStructure,
                                               bool persistent,
                                               void* comm
                                               )
    : AbstractGhostTimeCluster(maxTimeStepSize,
                               timeStepRate,
                               globalTimeClusterId,
                               otherGlobalTimeClusterId,
                               meshStructure), persistent(persistent),
                               comm(comm)
                               {
    if (persistent) {
      for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
        if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId) ) {
          MPI_Send_init(meshStructure->copyRegions[region],
                    static_cast<int>(meshStructure->copyRegionSizes[region]),
                    MPI_C_REAL,
                    meshStructure->neighboringClusters[region][0],
                    timeData + meshStructure->sendIdentifiers[region],
                    seissol::MPI::mpi.comm(),
                    meshStructure->sendRequests + region);
          MPI_Recv_init(meshStructure->ghostRegions[region],
                    static_cast<int>(meshStructure->ghostRegionSizes[region]),
                    MPI_C_REAL,
                    meshStructure->neighboringClusters[region][0],
                    timeData + meshStructure->receiveIdentifiers[region],
                    seissol::MPI::mpi.comm(),
                    meshStructure->receiveRequests + region);
        }
      }
    }
  }

  void DirectGhostTimeCluster::finalize() {
    AbstractGhostTimeCluster::finalize();
    if (persistent) {
      for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
        if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
          MPI_Request_free(meshStructure->sendRequests + region);
          MPI_Request_free(meshStructure->receiveRequests + region);
        }
      }
    }
  }
} // namespace seissol::time_stepping
