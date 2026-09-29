// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_CCLEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_CCLEXCHANGESCHEDULER_H_

#include "Solver/TimeStepping/Halo/Stream/StreamExchangeScheduler.h"

#include <cstddef>
#include <utility>
#include <vector>

namespace seissol::solver {

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

  /**
   * With `LaunchOrder::PerDirection`, sets up the connections of all communicators, one after the
   * other in the same order on all processes; collective over all processes.
   */
  void prepare() override;

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
  // for each slot: the transport that sends, and the one that receives in its direction
  std::vector<const ScheduledTransport*> senders_;
  std::vector<const ScheduledTransport*> receivers_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_CCLEXCHANGESCHEDULER_H_
