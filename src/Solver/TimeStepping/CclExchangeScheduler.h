// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLEXCHANGESCHEDULER_H_

#include "Solver/TimeStepping/ExchangeScheduler.h"

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace seissol::time_stepping {

/**
 * Exchanges the halo data with NCCL, RCCL or oneCCL, depending on the device backend. With
 * `LaunchOrder::Global`, all exchanges of a process go through one communicator and one stream;
 * with `LaunchOrder::PerDirection`, each direction between two time clusters that exchange data
 * has a communicator and a stream of its own.
 *
 * Only available in device builds with CCL support.
 */
class CclExchangeScheduler : public ExchangeScheduler {
  public:
  /**
   * Sets up the communicators; collective over all processes.
   */
  CclExchangeScheduler(std::size_t clusterCount, LaunchOrder order);
  ~CclExchangeScheduler() override;

  CclExchangeScheduler(const CclExchangeScheduler&) = delete;
  CclExchangeScheduler(CclExchangeScheduler&&) = delete;
  CclExchangeScheduler& operator=(const CclExchangeScheduler&) = delete;
  CclExchangeScheduler& operator=(CclExchangeScheduler&&) = delete;

  [[nodiscard]] void* latestEvent() const override { return latestEvent_; }
  void setLatestEvent(void* event) override { latestEvent_ = event; }
  [[nodiscard]] std::vector<void*> streams() const override;
  void releaseEvents() override;

  protected:
  void added(const ScheduledTransport& transport) override;
  Ticket launch(std::size_t from,
                std::size_t to,
                const ScheduledTransport* sender,
                const ScheduledTransport* receiver,
                const std::vector<void*>& after) override;
  bool completed(Ticket ticket) override;

  private:
  /// the communicator and stream of a direction
  [[nodiscard]] std::size_t index(std::size_t from, std::size_t to) const;

  std::size_t clusterCount_;
  std::vector<void*> communicators_;
  std::vector<void*> streams_;
  std::map<Ticket, void*> pendingEvents_;
  Ticket nextTicket_{0};

  // ordered on the device: the events of the groups launched since the last release, and the one of
  // the latest group
  std::vector<void*> launchedEvents_;
  void* latestEvent_{nullptr};
  std::vector<std::pair<void*, void*>> registrations_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLEXCHANGESCHEDULER_H_
