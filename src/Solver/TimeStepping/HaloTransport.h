// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORT_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORT_H_

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <memory>

namespace seissol::time_stepping {

/**
 * Moves the halo data between the copy layer of one time cluster and the ghost layer of one remote
 * time cluster: the copy regions go out, the ghost regions come in.
 *
 * Sending and receiving are independent of each other. Each of them is in flight at most once at a
 * time; it is started, and then tested until it has completed. When and how often that happens is
 * decided by the caller; the transport only moves the data.
 */
class HaloTransport {
  public:
  HaloTransport() = default;
  virtual ~HaloTransport() = default;

  HaloTransport(const HaloTransport&) = delete;
  HaloTransport(HaloTransport&&) = delete;
  HaloTransport& operator=(const HaloTransport&) = delete;
  HaloTransport& operator=(HaloTransport&&) = delete;

  /**
   * Starts sending the copy regions. They must not change until `testSend()` has returned true.
   */
  virtual void startSend() = 0;

  /**
   * Progresses the send. Returns true once it has completed, or if none is in flight.
   */
  virtual bool testSend() = 0;

  /**
   * Starts receiving into the ghost regions. They must not be read until `testReceive()` has
   * returned true.
   */
  virtual void startReceive() = 0;

  /**
   * Progresses the receive. Returns true once the data is in the ghost regions, or if none is in
   * flight.
   */
  virtual bool testReceive() = 0;

  /**
   * Releases everything that needs MPI; called before MPI is finalized.
   */
  virtual void finalize() {}
};

/**
 * Creates the transport for the given regions and transfer mode. Persistent transports set up their
 * MPI requests once and restart them for each transfer.
 */
std::unique_ptr<HaloTransport> createHaloTransport(const solver::RemoteClusterPair& regions,
                                                   Mpi::DataTransferMode mode,
                                                   bool persistent);

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORT_H_
