// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_SHMEMEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_SHMEMEXCHANGESCHEDULER_H_

#include "Solver/TimeStepping/StreamExchangeScheduler.h"

#include <cstddef>
#include <cstdint>
#include <map>
#include <tuple>
#include <vector>

namespace seissol::time_stepping {

/**
 * Exchanges the halo data with one-sided puts of NVSHMEM, ROCSHMEM or Intel SHMEM, depending on the
 * device backend, all in the global order on one stream.
 *
 * Each process has a staging window in symmetric memory, which holds the ghost regions of all its
 * transports; a sender puts its copy regions right there, at offsets the receivers have told it.
 * Per direction and peer, two signals count the exchanges: a receiver signals that its staging
 * window may be written ("clear to send"), a sender signals the arrival of its data. A group of an
 * exchange first signals clear to send to all peers it receives from, then waits for the clear to
 * send of the peers it sends to and puts its data, waits until the puts have left its buffers,
 * waits for the arrival of the data it receives, and copies it from the staging window into the
 * ghost regions. Since a group signals before it waits for anything, it can complete once all
 * peers have started their groups of the same exchange, as for two-sided operations.
 *
 * Only available in device builds with SHMEM support.
 */
class ShmemExchangeScheduler : public StreamExchangeScheduler {
  public:
  /**
   * Initializes the SHMEM library; collective over all processes.
   */
  explicit ShmemExchangeScheduler(std::size_t clusterCount);
  ~ShmemExchangeScheduler() override;

  ShmemExchangeScheduler(const ShmemExchangeScheduler&) = delete;
  ShmemExchangeScheduler(ShmemExchangeScheduler&&) = delete;
  ShmemExchangeScheduler& operator=(const ShmemExchangeScheduler&) = delete;
  ShmemExchangeScheduler& operator=(ShmemExchangeScheduler&&) = delete;

  /**
   * Allocates the staging window and the signals, and tells the senders where to put their data;
   * collective over all processes.
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
  [[nodiscard]] std::size_t signalIndex(std::size_t from, std::size_t to, int peer) const;

  std::vector<const ScheduledTransport*> transports_;
  // for each transport: where its ghost regions lie in the own staging window, and where its copy
  // regions go in the staging windows of the receivers
  std::map<const ScheduledTransport*, std::vector<std::uint64_t>> ghostOffsets_;
  std::map<const ScheduledTransport*, std::vector<std::uint64_t>> remoteOffsets_;
  std::size_t windowSize_{0};

  int rank_{0};
  int size_{1};
  char* window_{nullptr};
  std::uint64_t* clearToSend_{nullptr};
  std::uint64_t* arrived_{nullptr};
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_SHMEMEXCHANGESCHEDULER_H_
