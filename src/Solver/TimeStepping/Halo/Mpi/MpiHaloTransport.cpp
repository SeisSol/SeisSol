// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MpiHaloTransport.h"

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"

#include <cstddef>
#include <list>
#include <mpi.h>

namespace seissol::time_stepping {

bool testRequests(MPI_Request* requests, std::list<std::size_t>& regions) {
  for (auto region = regions.begin(); region != regions.end();) {
    MPI_Request* request = &requests[*region];
    int testSuccess = 0;
    MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess != 0) {
      region = regions.erase(region);
    } else {
      ++region;
    }
  }
  return regions.empty();
}

MpiHaloTransport::MpiHaloTransport(const solver::RemoteClusterPair& regions, bool persistent)
    : regions_(regions), persistent_(persistent), sendRequests_(regions.copy.size()),
      recvRequests_(regions.ghost.size()) {
  if (persistent_) {
    for (std::size_t region = 0; region < regions_.copy.size(); ++region) {
      MPI_Send_init(regions_.copy[region].data,
                    static_cast<int>(regions_.copy[region].size),
                    Mpi::precisionToMpiType(regions_.copy[region].datatype),
                    regions_.copy[region].rank,
                    regions_.copy[region].tag,
                    seissol::Mpi::mpi.comm(),
                    sendRequests_.data() + region);
    }
    for (std::size_t region = 0; region < regions_.ghost.size(); ++region) {
      MPI_Recv_init(regions_.ghost[region].data,
                    static_cast<int>(regions_.ghost[region].size),
                    Mpi::precisionToMpiType(regions_.ghost[region].datatype),
                    regions_.ghost[region].rank,
                    regions_.ghost[region].tag,
                    seissol::Mpi::mpi.comm(),
                    recvRequests_.data() + region);
    }
  }
}

void MpiHaloTransport::startSend() {
  if (persistent_) {
    MPI_Startall(static_cast<int>(sendRequests_.size()), sendRequests_.data());
  }
  for (std::size_t region = 0; region < regions_.copy.size(); ++region) {
    if (!persistent_) {
      MPI_Isend(regions_.copy[region].data,
                static_cast<int>(regions_.copy[region].size),
                Mpi::precisionToMpiType(regions_.copy[region].datatype),
                regions_.copy[region].rank,
                regions_.copy[region].tag,
                seissol::Mpi::mpi.comm(),
                sendRequests_.data() + region);
    }
    sendQueue_.push_back(region);
  }
}

bool MpiHaloTransport::testSend() { return testRequests(sendRequests_.data(), sendQueue_); }

void MpiHaloTransport::startReceive() {
  if (persistent_) {
    MPI_Startall(static_cast<int>(recvRequests_.size()), recvRequests_.data());
  }
  for (std::size_t region = 0; region < regions_.ghost.size(); ++region) {
    if (!persistent_) {
      MPI_Irecv(regions_.ghost[region].data,
                static_cast<int>(regions_.ghost[region].size),
                Mpi::precisionToMpiType(regions_.ghost[region].datatype),
                regions_.ghost[region].rank,
                regions_.ghost[region].tag,
                seissol::Mpi::mpi.comm(),
                recvRequests_.data() + region);
    }
    receiveQueue_.push_back(region);
  }
}

bool MpiHaloTransport::testReceive() { return testRequests(recvRequests_.data(), receiveQueue_); }

void MpiHaloTransport::finalize() {
  if (persistent_) {
    for (auto& request : sendRequests_) {
      MPI_Request_free(&request);
    }
    for (auto& request : recvRequests_) {
      MPI_Request_free(&request);
    }
  }
}

} // namespace seissol::time_stepping
