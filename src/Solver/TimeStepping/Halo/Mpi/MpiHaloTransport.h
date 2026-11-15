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
#include <string>
#include <vector>

namespace seissol::solver {

/**
 * Where the messages of the regions of a transport go in the trace (SEISSOL_TRACE): a track for
 * each region, named after its peer, with a span from the start of the exchange to the completion
 * of its message as the transport saw it.
 */
struct RegionTrace {
  const std::vector<RemoteCluster>* clusters{nullptr};
  bool send{false};
  /// the start of the exchange, in microseconds of the trace; negative while it does not record
  double begin{-1};
};

/**
 * Tests the requests of the given regions, and removes the completed ones from the list (and puts
 * them into the trace, if given). Returns true once the list is empty.
 */
bool testRequests(MPI_Request* requests,
                  std::list<std::size_t>& regions,
                  const RegionTrace* trace = nullptr);

/// the start of an exchange for RegionTrace::begin: now if the trace records, else negative
double regionTraceBegin();

/// puts the message of a region, completed now, into the trace (if it records)
void traceRegion(const RegionTrace& trace, std::size_t region);

/**
 * The messages of the regions in flight, with the rank, the tag and the size of each.
 */
std::string describePending(const RemoteClusterPair& regions,
                            const std::list<std::size_t>& sends,
                            const std::list<std::size_t>& receives);

/**
 * Sends and receives the halo regions directly with MPI, from and into the memory they reside in.
 */
class MpiHaloTransport : public HaloTransport {
  public:
  MpiHaloTransport(const RemoteClusterPair& regions, bool persistent);
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
  [[nodiscard]] std::string describe() const override;

  private:
  RemoteClusterPair regions_;
  bool persistent_;
  std::vector<MPI_Request> sendRequests_;
  std::vector<MPI_Request> recvRequests_;
  std::list<std::size_t> sendQueue_;
  std::list<std::size_t> receiveQueue_;
  double sendBegin_{-1};
  double receiveBegin_{-1};
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_MPI_MPIHALOTRANSPORT_H_
