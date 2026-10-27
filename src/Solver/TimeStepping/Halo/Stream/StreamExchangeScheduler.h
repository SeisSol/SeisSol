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
  StreamExchangeScheduler(std::size_t clusterCount, LaunchOrder order, std::size_t configCount);
  ~StreamExchangeScheduler() override;

  StreamExchangeScheduler(const StreamExchangeScheduler&) = delete;
  StreamExchangeScheduler(StreamExchangeScheduler&&) = delete;
  StreamExchangeScheduler& operator=(const StreamExchangeScheduler&) = delete;
  StreamExchangeScheduler& operator=(StreamExchangeScheduler&&) = delete;

  [[nodiscard]] std::vector<parallel::runtime::StreamRuntime*> runtimes() const override;
  [[nodiscard]] EventUsage eventUsage() const override;

  protected:
  /// the slot of a direction
  [[nodiscard]] std::size_t slot(std::size_t from, std::size_t to) const;

  /// the slots of all directions between neighboring time clusters (of any configurations), each
  /// once; the same on all processes
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
  // one for each used slot
  std::vector<std::unique_ptr<parallel::runtime::StreamRuntime>> runtimes_;
  // not ordered on the device: the events of the groups the host has not seen complete yet
  std::map<Ticket, parallel::runtime::Event> pendingEvents_;
#ifdef ACL_DEVICE
  Ticket nextTicket_{0};
#endif
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMEXCHANGESCHEDULER_H_
