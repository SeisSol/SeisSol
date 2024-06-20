#include "Solver/Clustering/Communication/DirectGhostTimeCluster.h"
#include "Parallel/MPI.h"
#include <Initializer/BasicTypedefs.hpp>
#include <Initializer/preProcessorMacros.hpp>
#include <Kernels/precision.hpp>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/AbstractGhostTimeCluster.h>
#include <cassert>
#include <mpi.h>

namespace seissol::time_stepping {
void DirectGhostTimeCluster::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.time.at(ComputeStep::Correct) > lastSendTime);
  lastSendTime = ct.time.at(ComputeStep::Correct);
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      if (persistent) {
        MPI_Start(meshStructure->sendRequests + region);
      } else {
        MPI_Isend(meshStructure->copyRegions[region],
                  static_cast<int>(meshStructure->copyRegionSizes[region]),
                  MPI_C_REAL,
                  meshStructure->neighboringClusters[region][0],
                  timeData + meshStructure->sendIdentifiers[region],
                  seissol::MPI::mpi.comm(),
                  meshStructure->sendRequests + region);
      }
      sendQueue.push_back(region);
    }
  }
}

void DirectGhostTimeCluster::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.time.at(ComputeStep::Predict) >= lastSendTime);
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      if (persistent) {
        MPI_Start(meshStructure->receiveRequests + region);
      } else {
        MPI_Irecv(meshStructure->ghostRegions[region],
                  static_cast<int>(meshStructure->ghostRegionSizes[region]),
                  MPI_C_REAL,
                  meshStructure->neighboringClusters[region][0],
                  timeData + meshStructure->receiveIdentifiers[region],
                  seissol::MPI::mpi.comm(),
                  meshStructure->receiveRequests + region);
      }
      receiveQueue.push_back(region);
    }
  }
}

bool DirectGhostTimeCluster::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testQueue(meshStructure->receiveRequests, receiveQueue);
}

DirectGhostTimeCluster::DirectGhostTimeCluster(double maxTimeStepSize,
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
  if (persistent) {
    for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
      if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
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
