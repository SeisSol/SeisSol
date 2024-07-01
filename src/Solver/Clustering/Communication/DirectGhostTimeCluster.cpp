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
  if (persistent) {
    MPI_Startall(sendRequests.size(), sendRequests.data());
    for (std::size_t i = 0; i < copyClusters.size(); ++i) {
      sendQueue.push_back(i);
    }
  } else {
    for (std::size_t i = 0; i < copyClusters.size(); ++i) {
      const auto& cluster = copyClusters.at(i);
      MPI_Isend(cluster.data,
                static_cast<int>(cluster.size),
                cluster.datatype,
                cluster.rank,
                cluster.tag,
                seissol::MPI::mpi.comm(),
                &sendRequests[i]);
      sendQueue.push_back(i);
    }
  }
}

void DirectGhostTimeCluster::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.time.at(ComputeStep::Predict) >= lastSendTime);
  if (persistent) {
    MPI_Startall(recvRequests.size(), recvRequests.data());
    for (std::size_t i = 0; i < ghostClusters.size(); ++i) {
      recvQueue.push_back(i);
    }
  } else {
    for (std::size_t i = 0; i < ghostClusters.size(); ++i) {
      const auto& cluster = ghostClusters.at(i);
      MPI_Irecv(cluster.data,
                static_cast<int>(cluster.size),
                cluster.datatype,
                cluster.rank,
                cluster.tag,
                seissol::MPI::mpi.comm(),
                &recvRequests[i]);
      recvQueue.push_back(i);
    }
  }
}

bool DirectGhostTimeCluster::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testQueue(recvRequests, recvQueue);
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
    for (std::size_t i = 0; i < copyClusters.size(); ++i) {
      {
        const auto& cluster = copyClusters.at(i);
        MPI_Send_init(cluster.data,
                      static_cast<int>(cluster.size),
                      cluster.datatype,
                      cluster.rank,
                      cluster.tag,
                      seissol::MPI::mpi.comm(),
                      &sendRequests[i]);
      }
      {
        const auto& cluster = ghostClusters.at(i);
        MPI_Recv_init(cluster.data,
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

void DirectGhostTimeCluster::finalize() {
  if (persistent) {
    for (auto& request : sendRequests) {
      MPI_Request_free(&request);
    }
    for (auto& request : recvRequests) {
      MPI_Request_free(&request);
    }
  }
}
} // namespace seissol::time_stepping
