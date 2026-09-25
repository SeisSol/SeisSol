// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMMPIEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMMPIEXCHANGESCHEDULER_H_

#include "Solver/TimeStepping/Halo/Stream/StreamExchangeScheduler.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::time_stepping {

/**
 * Exchanges the halo data with MPI operations that run in the order of a device stream: with MPICH
 * 4.1 and newer, the sends and receives of a group are enqueued on a stream communicator
 * (MPIX_Stream), followed by a wait for all of them; with HPE Cray MPICH, they are enqueued on a
 * queue that gets started and waited for on the stream (MPIX_Queue). One stream communicator or
 * queue per slot. The messages match by the tags of the regions.
 *
 * Only available in device builds with stream-aware MPI.
 */
class StreamMpiExchangeScheduler : public StreamExchangeScheduler {
  public:
  /**
   * Sets up the stream communicators or queues; collective over all processes.
   */
  StreamMpiExchangeScheduler(std::size_t clusterCount, LaunchOrder order);
  ~StreamMpiExchangeScheduler() override;

  StreamMpiExchangeScheduler(const StreamMpiExchangeScheduler&) = delete;
  StreamMpiExchangeScheduler(StreamMpiExchangeScheduler&&) = delete;
  StreamMpiExchangeScheduler& operator=(const StreamMpiExchangeScheduler&) = delete;
  StreamMpiExchangeScheduler& operator=(StreamMpiExchangeScheduler&&) = delete;

  protected:
  void enqueueGroup(std::size_t slot,
                    std::size_t from,
                    std::size_t to,
                    std::size_t exchange,
                    const ScheduledTransport* sender,
                    const ScheduledTransport* receiver) override;

  private:
  struct Slot;
  std::vector<std::unique_ptr<Slot>> slots_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_STREAMMPIEXCHANGESCHEDULER_H_
