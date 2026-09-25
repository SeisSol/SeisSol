// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_MPI_STAGEDMPIHALOTRANSPORT_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_MPI_STAGEDMPIHALOTRANSPORT_H_

#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <Device/device.h>
#include <cstddef>
#include <list>
#include <mpi.h>
#include <vector>

namespace seissol::time_stepping {

/**
 * Sends and receives the halo regions with MPI through pinned host buffers: the copy regions are
 * copied from the device before sending, the ghost regions to the device after receiving. For MPI
 * libraries that cannot work on device memory.
 */
class StagedMpiHaloTransport : public HaloTransport {
  public:
  StagedMpiHaloTransport(const solver::RemoteClusterPair& regions, bool persistent);
  ~StagedMpiHaloTransport() override;

  StagedMpiHaloTransport(const StagedMpiHaloTransport&) = delete;
  StagedMpiHaloTransport(StagedMpiHaloTransport&&) = delete;
  StagedMpiHaloTransport& operator=(const StagedMpiHaloTransport&) = delete;
  StagedMpiHaloTransport& operator=(StagedMpiHaloTransport&&) = delete;

  void startSend() override;
  bool testSend() override;
  void startReceive() override;
  bool testReceive() override;
  void finalize() override;

  private:
  enum class ReceiveState { RequiresMpiTesting, RequiresCopyTesting, Ready };

  solver::RemoteClusterPair regions_;
  bool persistent_;
  std::vector<MPI_Request> sendRequests_;
  std::vector<MPI_Request> recvRequests_;
  std::list<std::size_t> sendQueue_;
  std::list<std::size_t> receiveQueue_;
  std::vector<void*> hostCopyRegions_;
  std::vector<void*> hostGhostRegions_;
  std::vector<void*> copyStreams_;
  std::vector<void*> ghostStreams_;
  std::vector<ReceiveState> receiveStates_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_MPI_STAGEDMPIHALOTRANSPORT_H_
