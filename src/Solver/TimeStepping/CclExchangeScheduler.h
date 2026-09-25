// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLEXCHANGESCHEDULER_H_

#include "Solver/TimeStepping/StreamExchangeScheduler.h"

#include <cstddef>
#include <utility>
#include <vector>

namespace seissol::time_stepping {

/**
 * Exchanges the halo data with NCCL, RCCL or oneCCL, depending on the device backend: one
 * communicator per slot.
 *
 * Only available in device builds with CCL support.
 */
class CclExchangeScheduler : public StreamExchangeScheduler {
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

  protected:
  void added(const ScheduledTransport& transport) override;
  void enqueueGroup(std::size_t slot,
                    std::size_t from,
                    std::size_t to,
                    std::size_t exchange,
                    const ScheduledTransport* sender,
                    const ScheduledTransport* receiver) override;

  private:
  std::vector<void*> communicators_;
  std::vector<std::pair<void*, void*>> registrations_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLEXCHANGESCHEDULER_H_
