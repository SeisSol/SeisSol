// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifdef ACL_DEVICE

#include "StagedMpiHaloTransport.h"

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/MpiHaloTransport.h"

#include <Device/device.h>
#include <cstddef>
#include <list>
#include <mpi.h>

namespace seissol::time_stepping {

StagedMpiHaloTransport::StagedMpiHaloTransport(const solver::RemoteClusterPair& regions,
                                               bool persistent)
    : regions_(regions), persistent_(persistent), sendRequests_(regions.copy.size()),
      recvRequests_(regions.ghost.size()), hostCopyRegions_(regions.copy.size()),
      hostGhostRegions_(regions.ghost.size()), copyStreams_(regions.copy.size()),
      ghostStreams_(regions.ghost.size()),
      receiveStates_(regions.ghost.size(), ReceiveState::RequiresMpiTesting) {
  for (std::size_t region = 0; region < regions_.copy.size(); ++region) {
    copyStreams_[region] = device_.api->createStream();
    hostCopyRegions_[region] = device_.api->allocPinnedMem(
        regions_.copy[region].size * sizeOfRealType(regions_.copy[region].datatype));
    if (persistent_) {
      MPI_Send_init(hostCopyRegions_[region],
                    static_cast<int>(regions_.copy[region].size),
                    Mpi::precisionToMpiType(regions_.copy[region].datatype),
                    regions_.copy[region].rank,
                    regions_.copy[region].tag,
                    seissol::Mpi::mpi.comm(),
                    sendRequests_.data() + region);
    }
  }
  for (std::size_t region = 0; region < regions_.ghost.size(); ++region) {
    ghostStreams_[region] = device_.api->createStream();
    hostGhostRegions_[region] = device_.api->allocPinnedMem(
        regions_.ghost[region].size * sizeOfRealType(regions_.ghost[region].datatype));
    if (persistent_) {
      MPI_Recv_init(hostGhostRegions_[region],
                    static_cast<int>(regions_.ghost[region].size),
                    Mpi::precisionToMpiType(regions_.ghost[region].datatype),
                    regions_.ghost[region].rank,
                    regions_.ghost[region].tag,
                    seissol::Mpi::mpi.comm(),
                    recvRequests_.data() + region);
    }
  }
}

StagedMpiHaloTransport::~StagedMpiHaloTransport() {
  for (std::size_t region = 0; region < regions_.copy.size(); ++region) {
    device_.api->destroyGenericStream(copyStreams_[region]);
    device_.api->freePinnedMem(hostCopyRegions_[region]);
  }
  for (std::size_t region = 0; region < regions_.ghost.size(); ++region) {
    device_.api->destroyGenericStream(ghostStreams_[region]);
    device_.api->freePinnedMem(hostGhostRegions_[region]);
  }
}

void StagedMpiHaloTransport::startSend() {
  std::list<std::size_t> copying;
  for (std::size_t region = 0; region < regions_.copy.size(); ++region) {
    device_.api->copyFromAsync(hostCopyRegions_[region],
                               regions_.copy[region].data,
                               regions_.copy[region].size *
                                   sizeOfRealType(regions_.copy[region].datatype),
                               copyStreams_[region]);
    copying.push_back(region);
  }

  // send each region as soon as its data has arrived on the host
  while (!copying.empty()) {
    for (auto region = copying.begin(); region != copying.end();) {
      if (device_.api->isStreamWorkDone(copyStreams_[*region])) {
        if (persistent_) {
          MPI_Start(sendRequests_.data() + *region);
        } else {
          MPI_Isend(hostCopyRegions_[*region],
                    static_cast<int>(regions_.copy[*region].size),
                    Mpi::precisionToMpiType(regions_.copy[*region].datatype),
                    regions_.copy[*region].rank,
                    regions_.copy[*region].tag,
                    seissol::Mpi::mpi.comm(),
                    sendRequests_.data() + *region);
        }
        sendQueue_.push_back(*region);
        region = copying.erase(region);
      } else {
        ++region;
      }
    }
  }
}

bool StagedMpiHaloTransport::testSend() { return testRequests(sendRequests_.data(), sendQueue_); }

void StagedMpiHaloTransport::startReceive() {
  if (persistent_) {
    MPI_Startall(static_cast<int>(recvRequests_.size()), recvRequests_.data());
  }
  for (std::size_t region = 0; region < regions_.ghost.size(); ++region) {
    if (!persistent_) {
      MPI_Irecv(hostGhostRegions_[region],
                static_cast<int>(regions_.ghost[region].size),
                Mpi::precisionToMpiType(regions_.ghost[region].datatype),
                regions_.ghost[region].rank,
                regions_.ghost[region].tag,
                seissol::Mpi::mpi.comm(),
                recvRequests_.data() + region);
    }
    receiveStates_[region] = ReceiveState::RequiresMpiTesting;
    receiveQueue_.push_back(region);
  }
}

bool StagedMpiHaloTransport::testReceive() {
  for (auto region = receiveQueue_.begin(); region != receiveQueue_.end();) {
    switch (receiveStates_[*region]) {
    case ReceiveState::RequiresMpiTesting: {
      int testSuccess = 0;
      MPI_Test(recvRequests_.data() + *region, &testSuccess, MPI_STATUS_IGNORE);
      if (testSuccess != 0) {
        device_.api->copyToAsync(regions_.ghost[*region].data,
                                 hostGhostRegions_[*region],
                                 regions_.ghost[*region].size *
                                     sizeOfRealType(regions_.ghost[*region].datatype),
                                 ghostStreams_[*region]);
        receiveStates_[*region] = ReceiveState::RequiresCopyTesting;
      }
      ++region;
      break;
    }
    case ReceiveState::RequiresCopyTesting: {
      if (device_.api->isStreamWorkDone(ghostStreams_[*region])) {
        receiveStates_[*region] = ReceiveState::Ready;
        region = receiveQueue_.erase(region);
      } else {
        ++region;
      }
      break;
    }
    case ReceiveState::Ready: {
      region = receiveQueue_.erase(region);
      break;
    }
    }
  }
  return receiveQueue_.empty();
}

void StagedMpiHaloTransport::finalize() {
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

#endif // ACL_DEVICE
