// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_EXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_EXCHANGESCHEDULER_H_

#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/HaloTransport.h"

#include <cstddef>
#include <vector>

namespace seissol::time_stepping {

class ScheduledTransport;

/**
 * Launches the halo exchanges of one process for libraries whose point-to-point operations occupy
 * the stream they run on until the peers have posted their counterparts, as NCCL, RCCL and oneCCL
 * do.
 *
 * Each direction of exchange, from the copy layers of one time cluster to the ghost layers of
 * another one, has a communicator and a stream of its own. A process can take part in a direction
 * both as sender (with its copy layer of the first cluster) and as receiver (into its ghost layer
 * next to its copy layer of the second cluster). For each exchange, the scheduler waits until the
 * process is ready for all of its operations in that direction, and launches them as one group.
 * Consequently, all processes launch the groups of a direction in the same order, and a group
 * never waits for an operation that a process has queued behind another one.
 */
class ExchangeScheduler {
  public:
  using Ticket = std::size_t;

  explicit ExchangeScheduler(std::size_t clusterCount);
  virtual ~ExchangeScheduler() = default;

  ExchangeScheduler(const ExchangeScheduler&) = delete;
  ExchangeScheduler(ExchangeScheduler&&) = delete;
  ExchangeScheduler& operator=(const ExchangeScheduler&) = delete;
  ExchangeScheduler& operator=(ExchangeScheduler&&) = delete;

  /**
   * Makes the transport the sender from its cluster to its other cluster, and the receiver in the
   * opposite direction. Each direction has at most one sender and one receiver per process.
   */
  void add(ScheduledTransport& transport);

  /**
   * Marks the next send of the transport as ready; returns the index of its exchange.
   */
  std::size_t readySend(const ScheduledTransport& transport);

  /**
   * Marks the next receive of the transport as ready; returns the index of its exchange.
   */
  std::size_t readyReceive(const ScheduledTransport& transport);

  [[nodiscard]] bool sendCompleted(const ScheduledTransport& transport, std::size_t exchange);
  [[nodiscard]] bool receiveCompleted(const ScheduledTransport& transport, std::size_t exchange);

  protected:
  /**
   * Called once a transport has been added.
   */
  virtual void added(const ScheduledTransport& /*transport*/) {}

  /**
   * Launches one group of operations from cluster `from` to cluster `to`: the sends of `sender`
   * and the receives of `receiver`, either of which may be null.
   */
  virtual Ticket launch(std::size_t from,
                        std::size_t to,
                        const ScheduledTransport* sender,
                        const ScheduledTransport* receiver) = 0;

  /**
   * Whether the group behind the ticket has completed.
   */
  virtual bool completed(Ticket ticket) = 0;

  private:
  struct Direction {
    const ScheduledTransport* sender{nullptr};
    const ScheduledTransport* receiver{nullptr};
    std::size_t readySends{0};
    std::size_t readyReceives{0};
    std::vector<Ticket> groups;
  };

  Direction& direction(std::size_t from, std::size_t to);
  void launchReady(std::size_t from, std::size_t to);
  bool groupCompleted(Direction& direction, std::size_t exchange);

  std::size_t clusterCount_;
  std::vector<Direction> directions_;
};

/**
 * The transport between a copy layer and one ghost layer, as seen by an `ExchangeScheduler`: it
 * only marks its sends and receives as ready, the scheduler decides when they go out.
 */
class ScheduledTransport : public HaloTransport {
  public:
  ScheduledTransport(ExchangeScheduler& scheduler,
                     const solver::RemoteClusterPair& regions,
                     std::size_t cluster,
                     std::size_t otherCluster);

  void startSend() override;
  bool testSend() override;
  void startReceive() override;
  bool testReceive() override;

  [[nodiscard]] const solver::RemoteClusterPair& regions() const { return regions_; }
  [[nodiscard]] std::size_t cluster() const { return cluster_; }
  [[nodiscard]] std::size_t otherCluster() const { return otherCluster_; }

  private:
  ExchangeScheduler& scheduler_;
  solver::RemoteClusterPair regions_;
  std::size_t cluster_;
  std::size_t otherCluster_;
  bool sending_{false};
  bool receiving_{false};
  std::size_t sendExchange_{0};
  std::size_t receiveExchange_{0};
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_EXCHANGESCHEDULER_H_
