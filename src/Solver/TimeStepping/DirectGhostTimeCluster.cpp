// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Solver/TimeStepping/DirectGhostTimeCluster.h"

#include "Monitoring/Instrumentation.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"

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
  for (std::size_t region = 0; region < meshStructure.copy.size(); ++region) {
    if (!persistent) {
      MPI_Isend(meshStructure.copy[region].data,
                static_cast<int>(meshStructure.copy[region].size),
                MPI::precisionToMpiType(meshStructure.copy[region].datatype),
                meshStructure.copy[region].rank,
                meshStructure.copy[region].tag,
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
  for (std::size_t region = 0; region < meshStructure.ghost.size(); ++region) {
    if (!persistent) {
      MPI_Irecv(meshStructure.ghost[region].data,
                static_cast<int>(meshStructure.ghost[region].size),
                MPI::precisionToMpiType(meshStructure.ghost[region].datatype),
                meshStructure.ghost[region].rank,
                meshStructure.ghost[region].tag,
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
    for (std::size_t region = 0; region < this->meshStructure.copy.size(); ++region) {
      MPI_Send_init(this->meshStructure.copy[region].data,
                    static_cast<int>(this->meshStructure.copy[region].size),
                    MPI::precisionToMpiType(this->meshStructure.copy[region].datatype),
                    this->meshStructure.copy[region].rank,
                    this->meshStructure.copy[region].tag,
                    seissol::MPI::mpi.comm(),
                    sendRequests.data() + region);
    }
    for (std::size_t region = 0; region < this->meshStructure.ghost.size(); ++region) {
      MPI_Recv_init(this->meshStructure.ghost[region].data,
                    static_cast<int>(this->meshStructure.ghost[region].size),
                    MPI::precisionToMpiType(this->meshStructure.ghost[region].datatype),
                    this->meshStructure.ghost[region].rank,
                    this->meshStructure.ghost[region].tag,
                    seissol::MPI::mpi.comm(),
                    recvRequests.data() + region);
    }
  }
}

void DirectGhostTimeCluster::finalize() {
  if (persistent) {
    for (std::size_t region = 0; region < this->meshStructure.copy.size(); ++region) {
      MPI_Request_free(sendRequests.data() + region);
    }
    for (std::size_t region = 0; region < this->meshStructure.ghost.size(); ++region) {
      MPI_Request_free(recvRequests.data() + region);
    }
  }
}
} // namespace seissol::time_stepping
