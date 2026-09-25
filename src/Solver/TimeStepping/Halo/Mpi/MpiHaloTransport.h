// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_MPI_MPIHALOTRANSPORT_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_MPI_MPIHALOTRANSPORT_H_

#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <cstddef>
#include <list>
#include <mpi.h>
#include <vector>

namespace seissol::time_stepping {

/**
 * Tests the requests of the given regions, and removes the completed ones from the list. Returns
 * true once the list is empty.
 */
bool testRequests(MPI_Request* requests, std::list<std::size_t>& regions);

/**
 * Sends and receives the halo regions directly with MPI, from and into the memory they reside in.
 */
class MpiHaloTransport : public HaloTransport {
  public:
  MpiHaloTransport(const solver::RemoteClusterPair& regions, bool persistent);
  ~MpiHaloTransport() override = default;

  MpiHaloTransport(const MpiHaloTransport&) = delete;
  MpiHaloTransport(MpiHaloTransport&&) = delete;
  MpiHaloTransport& operator=(const MpiHaloTransport&) = delete;
  MpiHaloTransport& operator=(MpiHaloTransport&&) = delete;

  void startSend() override;
  bool testSend() override;
  void startReceive() override;
  bool testReceive() override;
  void finalize() override;

  private:
  solver::RemoteClusterPair regions_;
  bool persistent_;
  std::vector<MPI_Request> sendRequests_;
  std::vector<MPI_Request> recvRequests_;
  std::list<std::size_t> sendQueue_;
  std::list<std::size_t> receiveQueue_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_MPI_MPIHALOTRANSPORT_H_
