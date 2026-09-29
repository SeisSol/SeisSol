// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMEXCHANGESCHEDULER_H_

#include "Parallel/Runtime/EventPool.h"
#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"

#include <cstddef>
#include <map>
#include <vector>

namespace seissol::solver {

/**
 * An exchange scheduler whose library enqueues the operations of a group on a device stream (CCL,
 * stream-aware MPI, SHMEM). With `LaunchOrder::Global`, all directions share one slot (and thus one
 * stream); otherwise, each direction between neighboring time clusters has a slot of its own.
 *
 * Before a group, its stream waits for the events the group has to wait for; after it, an event
 * gets recorded. The host either tests that event, or, ordered on the device, the dependent work
 * waits for it. The events come from a pool that hands out an event again only once nobody refers
 * to it any more: the scheduler until it has seen the group complete, or until a later group is the
 * latest one; and everybody who has got it as latestEvent() until they have waited for it.
 *
 * Only available in device builds.
 */
class StreamExchangeScheduler : public ExchangeScheduler {
  public:
  StreamExchangeScheduler(std::size_t clusterCount, LaunchOrder order);
  ~StreamExchangeScheduler() override;

  StreamExchangeScheduler(const StreamExchangeScheduler&) = delete;
  StreamExchangeScheduler(StreamExchangeScheduler&&) = delete;
  StreamExchangeScheduler& operator=(const StreamExchangeScheduler&) = delete;
  StreamExchangeScheduler& operator=(StreamExchangeScheduler&&) = delete;

  [[nodiscard]] std::vector<void*> streams() const override;
  [[nodiscard]] EventUsage eventUsage() const override;

  protected:
  /// the slot of a direction
  [[nodiscard]] std::size_t slot(std::size_t from, std::size_t to) const;

  /// the slots of all directions between neighboring time clusters, each once
  [[nodiscard]] std::vector<std::size_t> usedSlots() const;

  [[nodiscard]] void* stream(std::size_t slot) const { return streams_[slot]; }

  [[nodiscard]] std::size_t clusterCount() const { return clusterCount_; }

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

  private:
  Ticket launch(std::size_t from,
                std::size_t to,
                std::size_t exchange,
                const ScheduledTransport* sender,
                const ScheduledTransport* receiver,
                const std::vector<void*>& after) final;
  bool completed(Ticket ticket) final;

  std::size_t clusterCount_;
  std::vector<void*> streams_;
  parallel::runtime::EventPool events_;
  // not ordered on the device: the events of the groups the host has not seen complete yet
  std::map<Ticket, parallel::runtime::EventRef> pendingEvents_;
  Ticket nextTicket_{0};
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMEXCHANGESCHEDULER_H_
