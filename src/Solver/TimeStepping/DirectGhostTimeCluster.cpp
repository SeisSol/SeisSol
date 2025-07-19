// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Solver/TimeStepping/DirectGhostTimeCluster.h"
#include "Parallel/MPI.h"
#include <Monitoring/Instrumentation.h>
#include <Solver/TimeStepping/AbstractGhostTimeCluster.h>
#include <Solver/TimeStepping/HaloCommunication.h>
#include <cassert>
#include <cstddef>
#include <mpi.h>

namespace seissol::time_stepping {
void DirectGhostTimeCluster::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.correctionTime > lastSendTime);
  lastSendTime = ct.correctionTime;
  if (persistent) {
    MPI_Startall(sendRequests.size(), sendRequests.data());
  }
  for (std::size_t region = 0; region < meshStructure.size(); ++region) {
    if (!persistent) {
      MPI_Isend(meshStructure[region].copy.data,
                static_cast<int>(meshStructure[region].copy.size),
                MPI::precisionToMpiType(meshStructure[region].copy.datatype),
                meshStructure[region].copy.rank,
                meshStructure[region].copy.tag,
                seissol::MPI::mpi.comm(),
                sendRequests.data() + region);
    }
    sendQueue.push_back(region);
  }
}

void DirectGhostTimeCluster::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.predictionTime >= lastSendTime);
  if (persistent) {
    MPI_Startall(recvRequests.size(), recvRequests.data());
  }
  for (std::size_t region = 0; region < meshStructure.size(); ++region) {
    if (!persistent) {
      MPI_Irecv(meshStructure[region].ghost.data,
                static_cast<int>(meshStructure[region].ghost.size),
                MPI::precisionToMpiType(meshStructure[region].ghost.datatype),
                meshStructure[region].ghost.rank,
                meshStructure[region].ghost.tag,
                seissol::MPI::mpi.comm(),
                recvRequests.data() + region);
    }
    receiveQueue.push_back(region);
  }
}

bool DirectGhostTimeCluster::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testQueue(recvRequests.data(), receiveQueue);
}

DirectGhostTimeCluster::DirectGhostTimeCluster(
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
  if (persistent) {
    for (std::size_t region = 0; region < this->meshStructure.size(); ++region) {
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

void DirectGhostTimeCluster::finalize() {
  if (persistent) {
    for (std::size_t region = 0; region < this->meshStructure.size(); ++region) {
      MPI_Request_free(sendRequests.data() + region);
      MPI_Request_free(recvRequests.data() + region);
    }
  }
}
} // namespace seissol::time_stepping
