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
#include <deque>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

namespace seissol::time_stepping {

class ScheduledTransport;

/**
 * In which order the groups of the different directions go out.
 */
enum class LaunchOrder {
  /// each direction on a stream of its own, independently of the others
  PerDirection,
  /// all directions on one stream, in an order that is the same on all processes
  Global
};

/**
 * Launches the halo exchanges of one process for libraries whose point-to-point operations occupy
 * the stream they run on until the peers have posted their counterparts, as NCCL, RCCL and oneCCL
 * do.
 *
 * A direction of exchange goes from the copy layers of one time cluster to the ghost layers of
 * another one. A process can take part in a direction both as sender (with its copy layer of the
 * first cluster) and as receiver (into its ghost layer next to its copy layer of the second
 * cluster). For each exchange, the scheduler waits until the process is ready for all of its
 * operations in that direction, and launches them as one group. Consequently, all processes
 * launch the groups of a direction in the same order.
 *
 * With `LaunchOrder::PerDirection`, each direction needs a stream of its own. With
 * `LaunchOrder::Global`, the groups of all directions go out on one stream, ordered by the point
 * in logical time at which their data is complete: the start of the prediction of the sending
 * cluster that completes it, counted in steps of the smallest cluster since the synchronization
 * point. Ties go by the clusters of the direction. Every process thus launches its groups in the
 * same order, and every step waits only for groups that come before the ones its own data goes
 * out with.
 */
class ExchangeScheduler {
  public:
  using Ticket = std::size_t;

  ExchangeScheduler(std::size_t clusterCount, LaunchOrder order);
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
   * Announces the exchanges of the transport up to the next synchronization point. With
   * `LaunchOrder::Global`, nothing goes out until all transports have announced them.
   */
  void startInterval(const ScheduledTransport& transport, const ExchangeInterval& interval);

  /**
   * Marks the next send of the transport as ready, once the work behind the event has completed
   * on the device; returns the index of its exchange.
   */
  std::size_t readySend(const ScheduledTransport& transport, void* after = nullptr);

  /**
   * Marks the next receive of the transport as ready, once the work behind the event has
   * completed on the device; returns the index of its exchange.
   */
  std::size_t readyReceive(const ScheduledTransport& transport, void* after = nullptr);

  /**
   * Orders the groups on the device: a group starts after the events its operations were made
   * ready with, and counts as done for the host as soon as it has been launched; the work that
   * depends on it waits for latestEvent(). Only with `LaunchOrder::Global`, for clusters that wait
   * for each other on the device.
   */
  void setStreamOrdered(bool streamOrdered);
  [[nodiscard]] bool streamOrdered() const { return streamOrdered_; }

  /**
   * Without launching, the groups only count as launched; their operations come from elsewhere,
   * e.g. a recording.
   */
  void setLaunching(bool launching) { launching_ = launching; }

  /**
   * Launches only the groups whose data is complete before the given point in logical time of the
   * current interval (with `LaunchOrder::Global`); the others wait until the horizon moves on.
   * Launches everything once it is back at its default.
   */
  void setHorizon(long time);

  /**
   * Whether all groups whose data is complete before the given point in logical time of the
   * current interval have been launched (with `LaunchOrder::Global`).
   */
  [[nodiscard]] bool launchedBefore(long time) const;

  /**
   * The event that completes with all groups launched so far, when ordered on the device; null if
   * there is none.
   */
  [[nodiscard]] virtual void* latestEvent() const { return nullptr; }

  /**
   * Sets up what needs all transports; collective over all processes, after all transports have
   * been added.
   */
  virtual void prepare() {}

  /**
   * Makes the given event stand for all groups launched so far, e.g. once they have run as part of
   * a replayed recording.
   */
  virtual void setLatestEvent(void* /*event*/) {}

  /**
   * The streams the groups run on.
   */
  [[nodiscard]] virtual std::vector<void*> streams() const { return {}; }

  /**
   * Releases what the groups launched so far needed; only once the device has completed them.
   */
  virtual void releaseEvents() {}

  [[nodiscard]] bool sendCompleted(const ScheduledTransport& transport, std::size_t exchange);
  [[nodiscard]] bool receiveCompleted(const ScheduledTransport& transport, std::size_t exchange);

  [[nodiscard]] LaunchOrder launchOrder() const { return order_; }

  protected:
  /**
   * Called once a transport has been added.
   */
  virtual void added(const ScheduledTransport& /*transport*/) {}

  /**
   * Launches the group of operations of the given exchange (counted per direction) from cluster
   * `from` to cluster `to`: the sends of `sender` and the receives of `receiver`, either of which
   * may be null, once the work behind the events `after` has completed.
   */
  virtual Ticket launch(std::size_t from,
                        std::size_t to,
                        std::size_t exchange,
                        const ScheduledTransport* sender,
                        const ScheduledTransport* receiver,
                        const std::vector<void*>& after) = 0;

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

    // the events each ready operation waits for, in the order of the exchanges
    std::deque<void*> sendsAfter;
    std::deque<void*> receivesAfter;

    // the sending cluster of the current interval
    long sendRate{1};
    long sendSteps{0};
    long exchangePeriod{1};
  };

  Direction& direction(std::size_t from, std::size_t to);
  [[nodiscard]] static bool ready(const Direction& direction);
  void launchReady(std::size_t from, std::size_t to);
  void launchNext(std::size_t from, std::size_t to);
  void launchInOrder();
  void orderInterval();
  bool groupCompleted(Direction& direction, std::size_t exchange);

  std::size_t clusterCount_;
  LaunchOrder order_;
  std::vector<Direction> directions_;
  std::size_t transports_{0};

  // for LaunchOrder::Global: the groups of the current interval (time, from, to), in order
  std::size_t announced_{0};
  std::vector<std::tuple<long, std::size_t, std::size_t>> sequence_;
  std::size_t launched_{0};

  bool streamOrdered_{false};
  bool launching_{true};
  long horizon_{std::numeric_limits<long>::max()};
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

  void startInterval(const ExchangeInterval& interval) override;
  [[nodiscard]] bool streamOrdered() const override { return scheduler_.streamOrdered(); }
  [[nodiscard]] void* latestEvent() const override { return scheduler_.latestEvent(); }
  void startSendAfter(void* event) override;
  void startReceiveAfter(void* event) override;
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
