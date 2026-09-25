// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_HALOTRANSPORT_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_HALOTRANSPORT_H_

#include "Solver/TimeStepping/Halo/HaloCommunication.h"

namespace seissol::solver {

/**
 * The exchanges of a transport up to the next synchronization point. Each direction is described
 * by its sending cluster: its time step rate, and its number of steps (of the smallest cluster) up
 * to the synchronization point. Both directions exchange once per exchange period.
 */
struct ExchangeInterval {
  long sendRate{1};
  long sendSteps{0};
  long receiveRate{1};
  long receiveSteps{0};
  long exchangePeriod{1};
};

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
   * Whether the transport orders its operations on the device: they start after the events they
   * were started with, and the work that depends on them waits for latestEvent(). testSend() and
   * testReceive() then report whether they have been launched, not whether they have completed.
   */
  [[nodiscard]] virtual bool streamOrdered() const { return false; }

  /**
   * For stream-ordered transports: an event that completes with all operations launched so far;
   * null if there is none.
   */
  [[nodiscard]] virtual void* latestEvent() const { return nullptr; }

  /**
   * Starts sending once the work behind the event has completed on the device (for stream-ordered
   * transports; the others start right away).
   */
  virtual void startSendAfter(void* /*event*/) { startSend(); }

  /**
   * Starts receiving once the work behind the event has completed on the device (for
   * stream-ordered transports; the others start right away).
   */
  virtual void startReceiveAfter(void* /*event*/) { startReceive(); }

  /**
   * Announces the exchanges up to the next synchronization point, before the first of them.
   */
  virtual void startInterval(const ExchangeInterval& /*interval*/) {}

  /**
   * Releases everything that needs MPI; called before MPI is finalized.
   */
  virtual void finalize() {}
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_HALOTRANSPORT_H_
