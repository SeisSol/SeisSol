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
  lastSendTime_ = ct_.correctionTime;
  if (persistent_) {
    MPI_Startall(sendRequests_.size(), sendRequests_.data());
  }
  for (std::size_t region = 0; region < meshStructure_.copy.size(); ++region) {
    if (!persistent_) {
      MPI_Isend(meshStructure_.copy[region].data,
                static_cast<int>(meshStructure_.copy[region].size),
                Mpi::precisionToMpiType(meshStructure_.copy[region].datatype),
                meshStructure_.copy[region].rank,
                meshStructure_.copy[region].tag,
                seissol::Mpi::mpi.comm(),
                sendRequests_.data() + region);
    }
    sendQueue_.push_back(region);
  }
}

void DirectGhostTimeCluster::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct.predictionTime >= lastSendTime);
  if (persistent_) {
    MPI_Startall(recvRequests_.size(), recvRequests_.data());
  }
  for (std::size_t region = 0; region < meshStructure_.ghost.size(); ++region) {
    if (!persistent_) {
      MPI_Irecv(meshStructure_.ghost[region].data,
                static_cast<int>(meshStructure_.ghost[region].size),
                Mpi::precisionToMpiType(meshStructure_.ghost[region].datatype),
                meshStructure_.ghost[region].rank,
                meshStructure_.ghost[region].tag,
                seissol::Mpi::mpi.comm(),
                recvRequests_.data() + region);
    }
    receiveQueue_.push_back(region);
  }
}

bool DirectGhostTimeCluster::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testQueue(recvRequests_.data(), receiveQueue_);
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
      persistent_(persistent) {
  if (persistent) {
    for (std::size_t region = 0; region < this->meshStructure_.copy.size(); ++region) {
      MPI_Send_init(this->meshStructure_.copy[region].data,
                    static_cast<int>(this->meshStructure_.copy[region].size),
                    Mpi::precisionToMpiType(this->meshStructure_.copy[region].datatype),
                    this->meshStructure_.copy[region].rank,
                    this->meshStructure_.copy[region].tag,
                    seissol::Mpi::mpi.comm(),
                    sendRequests_.data() + region);
    }
    for (std::size_t region = 0; region < this->meshStructure_.ghost.size(); ++region) {
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

void DirectGhostTimeCluster::finalize() {
  if (persistent_) {
    for (std::size_t region = 0; region < this->meshStructure_.copy.size(); ++region) {
      MPI_Request_free(sendRequests_.data() + region);
    }
    for (std::size_t region = 0; region < this->meshStructure_.ghost.size(); ++region) {
      MPI_Request_free(recvRequests_.data() + region);
    }
  }
}
} // namespace seissol::time_stepping
