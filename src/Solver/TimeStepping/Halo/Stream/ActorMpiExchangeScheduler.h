// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_ACTORMPIEXCHANGESCHEDULER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_ACTORMPIEXCHANGESCHEDULER_H_

#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/StreamExchangeScheduler.h"

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>

namespace seissol::solver {

/**
 * Exchanges the halo data with MPI in the order of a device stream, as stream-mpi does, but with
 * any MPI: from the host, through buffers in pinned host memory, by the progress thread of the host
 * actor (see parallel::runtime::HostActor::startProgressThread()).
 *
 * In a group, the stream first copies the copy regions of the sends into the host buffers of their
 * transport, and counts that in host memory; the thread then sends them. Then the stream waits
 * until the thread has received the ghost regions of the receives (into host buffers as well, once
 * the stream has copied the ones before out of them), and copies them into the ghost regions. The
 * thread counts the completed sends and receives in host memory, which the stream waits for (see
 * device::AbstractAPI::streamWaitMemory()); the copies are kernels (see
 * kernels::device::aux::hostcopy).
 *
 * So the host enqueues the work of the device without waiting for MPI; and neither the thread nor
 * MPI ever calls the device runtime or works on device memory: the thread that enqueues the work of
 * the device may block in the device runtime meanwhile, and hold its locks (e.g. in a launch, while
 * the device waits). A stream that waits for the host blocks the other streams on its hardware
 * queue as well (CUDA has 8 of them by default, see CUDA_DEVICE_MAX_CONNECTIONS): in the global
 * order, only the stream of the exchange waits (or its lanes). With CUDA, the kernels have to be
 * loaded before (CUDA_MODULE_LOADING=EAGER): loading one on its first launch may wait for the
 * device.
 *
 * Only available in device builds whose streams wait for host memory (CUDA, HIP).
 */
class ActorMpiExchangeScheduler : public StreamExchangeScheduler {
  public:
  /// with persistent MPI requests, if `persistent`
  ActorMpiExchangeScheduler(std::size_t clusterCount,
                            LaunchOrder order,
                            std::size_t configCount,
                            std::size_t lanes,
                            bool persistent);
  ~ActorMpiExchangeScheduler() override;

  ActorMpiExchangeScheduler(const ActorMpiExchangeScheduler&) = delete;
  ActorMpiExchangeScheduler(ActorMpiExchangeScheduler&&) = delete;
  ActorMpiExchangeScheduler& operator=(const ActorMpiExchangeScheduler&) = delete;
  ActorMpiExchangeScheduler& operator=(ActorMpiExchangeScheduler&&) = delete;

  /**
   * No: the host takes up the MPI operations of a group when it enqueues the group, not when a
   * replay runs it.
   */
  [[nodiscard]] bool recordable() const override { return false; }

  protected:
  void added(const ScheduledTransport& transport) override;
  void enqueueGroup(std::size_t slot,
                    std::size_t from,
                    std::size_t to,
                    std::size_t exchange,
                    const ScheduledTransport* sender,
                    const ScheduledTransport* receiver) override;

  private:
  struct Staging;
#ifdef ACL_DEVICE
  /// lets the stream wait until the counter has reached the value
  void waitFor(void* stream, std::uint32_t* counter, std::uint32_t value);

  /// enqueues a group without staging: MPI works on the regions on the device
  void enqueueUnstaged(void* deviceStream,
                       const ScheduledTransport* sender,
                       const ScheduledTransport* receiver);

  bool persistent_;
  // with a kernel of its own (see kernels::device::aux::hostcopy::waitHost()), instead of the wait
  // of the device
  bool waitKernel_{false};
  // through host buffers (by default), or MPI on the regions on the device
  bool staged_{true};
  // for each transport: its host buffers, its counters, and its MPI operations
  std::map<const ScheduledTransport*, std::unique_ptr<Staging>> staging_;
#endif
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_ACTORMPIEXCHANGESCHEDULER_H_
