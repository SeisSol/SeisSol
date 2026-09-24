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

#include <cstddef>
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

class ExchangeScheduler;

/**
 * Creates the transports of one process, and holds what they share.
 */
class HaloTransportFactory {
  public:
  /**
   * Persistent MPI transports set up their requests once and restart them for each transfer. The
   * CCL transports set up their communicators here, collectively over all processes.
   */
  HaloTransportFactory(Mpi::DataTransferMode mode, bool persistent, std::size_t clusterCount);
  ~HaloTransportFactory();

  HaloTransportFactory(const HaloTransportFactory&) = delete;
  HaloTransportFactory(HaloTransportFactory&&) = delete;
  HaloTransportFactory& operator=(const HaloTransportFactory&) = delete;
  HaloTransportFactory& operator=(HaloTransportFactory&&) = delete;

  /**
   * Creates the transport between the copy layer of `cluster` and the ghost layer of
   * `otherCluster`. The transports need to be destroyed before the factory.
   */
  std::unique_ptr<HaloTransport> create(const solver::RemoteClusterPair& regions,
                                        std::size_t cluster,
                                        std::size_t otherCluster);

  private:
  Mpi::DataTransferMode mode_;
  bool persistent_;
  std::unique_ptr<ExchangeScheduler> scheduler_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORT_H_
