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
 * Exchanges the halo data with NCCL, RCCL or oneCCL, depending on the device backend. Each
 * direction between two time clusters that exchange data has a communicator and a stream of its
 * own.
 *
 * Only available in device builds with CCL support.
 */
class CclExchangeScheduler : public ExchangeScheduler {
  public:
  /**
   * Sets up the communicators; collective over all processes.
   */
  explicit CclExchangeScheduler(std::size_t clusterCount);
  ~CclExchangeScheduler() override;

  CclExchangeScheduler(const CclExchangeScheduler&) = delete;
  CclExchangeScheduler(CclExchangeScheduler&&) = delete;
  CclExchangeScheduler& operator=(const CclExchangeScheduler&) = delete;
  CclExchangeScheduler& operator=(CclExchangeScheduler&&) = delete;

  protected:
  void added(const ScheduledTransport& transport) override;
  Ticket launch(std::size_t from,
                std::size_t to,
                const ScheduledTransport* sender,
                const ScheduledTransport* receiver) override;
  bool completed(Ticket ticket) override;

  private:
  [[nodiscard]] std::size_t index(std::size_t from, std::size_t to) const;

  std::size_t clusterCount_;
  std::vector<void*> communicators_;
  std::vector<void*> streams_;
  std::map<Ticket, void*> pendingEvents_;
  Ticket nextTicket_{0};
  std::vector<std::pair<void*, void*>> registrations_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLEXCHANGESCHEDULER_H_
