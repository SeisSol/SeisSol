// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMEXCHANGESCHEDULER_H_

#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"

#include <cstddef>
#include <map>
#include <memory>
#include <vector>

namespace seissol::solver {

/**
 * An exchange scheduler whose library enqueues the operations of a group on a device stream (CCL,
 * stream-aware MPI, SHMEM). With `LaunchOrder::Global`, all directions share one slot (and thus one
 * stream runtime); otherwise, each direction between neighboring time clusters has a slot of its
 * own.
 *
 * In the global order, the groups may go to several lanes instead (experimental): slots of streams
 * of the highest priority (which HIP puts on hardware queues of their own, apart from the work of
 * the clusters), where the groups go one after another as before, but the lanes run at the same
 * time. A group goes to a lane by its direction and its exchange, the same on all processes; each
 * lane thus gets its groups in the global order everywhere, so that operations that occupy their
 * stream until the peer has posted its counterpart cannot wait for each other in a cycle across
 * lanes either. Ordered on the device, the latest event then joins all lanes, on a stream of its
 * own.
 *
 * Before a group, its stream waits for the events the group has to wait for; after it, an event
 * gets recorded. The host either tests that event, or, ordered on the device, the dependent work
 * waits for it. The events come from the pool of the runtime, which hands out an event again only
 * once nobody refers to it any more: the scheduler until it has seen the group complete, or until a
 * later group is the latest one; and everybody who has got it as latestEvent() until they have
 * waited for it. While a super-timestep gets recorded, the groups become nodes of the recording
 * instead, and their events these nodes (see parallel::runtime::StreamRuntime::attach()).
 *
 * Only available in device builds.
 */
class StreamExchangeScheduler : public ExchangeScheduler {
  public:
  /// with `lanes` lanes in the global order (see above), one by default; `libraryStreams`: the
  /// library needs streams of the device (CCL, stream-aware MPI, SHMEM), which then stay streams of
  /// HIP with own queues as well (see parallel::runtime::NoOwnQueue)
  StreamExchangeScheduler(std::size_t clusterCount,
                          LaunchOrder order,
                          std::size_t configCount,
                          std::size_t lanes = 1,
                          bool libraryStreams = true);
  ~StreamExchangeScheduler() override;

  StreamExchangeScheduler(const StreamExchangeScheduler&) = delete;
  StreamExchangeScheduler(StreamExchangeScheduler&&) = delete;
  StreamExchangeScheduler& operator=(const StreamExchangeScheduler&) = delete;
  StreamExchangeScheduler& operator=(StreamExchangeScheduler&&) = delete;

  [[nodiscard]] std::vector<parallel::runtime::StreamRuntime*> runtimes() const override;
  [[nodiscard]] EventUsage eventUsage() const override;

  protected:
  /// the slot of a direction; the first one with lanes
  [[nodiscard]] std::size_t slot(std::size_t from, std::size_t to) const;

  /// the slot of exchange `exchange` of a direction: that of the direction, or its lane
  [[nodiscard]] std::size_t slot(std::size_t from, std::size_t to, std::size_t exchange) const;

  /// the number of slots (used or not), which the slots index
  [[nodiscard]] std::size_t slotCount() const;

  /// the slots of all directions between neighboring time clusters (of any configurations), each
  /// once, or all lanes; the same on all processes
  [[nodiscard]] std::vector<std::size_t> usedSlots() const;

  /// the stream to enqueue the operations of the slot on (see
  /// parallel::runtime::StreamRuntime::stream())
  [[nodiscard]] void* stream(std::size_t slot) { return runtimes_[slot]->stream(); }

  [[nodiscard]] std::size_t clusterCount() const { return clusterCount_; }

  /**
   * Captures the groups into graphs of their own when super-timesteps get recorded (see
   * parallel::runtime::StreamRuntime::isolateNodes()).
   */
  void isolateGroups();

  /**
   * Enqueues the operations of a group on the stream of its slot.
   */
  virtual void enqueueGroup(std::size_t slot,
                            std::size_t from,
                            std::size_t to,
                            std::size_t exchange,
                            const ScheduledTransport* sender,
                            const ScheduledTransport* receiver) = 0;

  /**
   * Waits until the device has completed all groups; before the library releases its resources.
   */
  void synchronize();

  Ticket launch(std::size_t from,
                std::size_t to,
                std::size_t exchange,
                const ScheduledTransport* sender,
                const ScheduledTransport* receiver,
                const std::vector<ActorEvent>& after) final;
  bool completed(Ticket ticket) final;

  private:
  std::size_t clusterCount_;
  std::size_t lanes_;
  // one for each used slot
  std::vector<std::unique_ptr<parallel::runtime::StreamRuntime>> runtimes_;
  // with several lanes: where the latest event joins them
  std::unique_ptr<parallel::runtime::StreamRuntime> join_;
  // not ordered on the device: the events of the groups the host has not seen complete yet
  std::map<Ticket, parallel::runtime::Event> pendingEvents_;
#ifdef ACL_DEVICE
  Ticket nextTicket_{0};
#endif
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMEXCHANGESCHEDULER_H_
