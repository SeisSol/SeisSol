// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_TASKHALOTRANSPORT_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_TASKHALOTRANSPORT_H_

#include "Parallel/Runtime/HostRuntime.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <atomic>
#include <cstddef>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#ifdef SEISSOL_HOST_TASK_GRAPH

#include <omp.h>

namespace seissol::solver {

/**
 * Starts and tests the exchanges of the halo exchange in a task graph (see TaskHaloTransport), and
 * completes their tasks once the exchanges have completed. The tasks hand their exchanges over from
 * any thread; one thread at a time polls.
 */
class TaskHaloTransport;

class TaskHaloEngine {
  public:
  /**
   * Counts a task of an exchange; before the task gets created.
   */
  void created() { outstanding_.fetch_add(1, std::memory_order_relaxed); }

  /**
   * From the task of an exchange of `owner`, once the work before it is done: the exchange may
   * start; the task completes with `event`.
   */
  void enqueue(const TaskHaloTransport* owner,
               HaloTransport* transport,
               bool send,
               omp_event_handle_t event);

  /**
   * Starts the exchanges handed over, in that order, tests those in flight, and completes the tasks
   * of those that have completed.
   */
  void poll();

  /// the tasks of exchanges that have been created and have not completed yet
  [[nodiscard]] std::size_t outstanding() const {
    return outstanding_.load(std::memory_order_acquire);
  }

  [[nodiscard]] std::string describe() const;

  private:
  struct Exchange {
    const TaskHaloTransport* owner{nullptr};
    HaloTransport* transport{nullptr};
    bool send{false};
    // NOLINTNEXTLINE(bugprone-invalid-enum-default-initialization)
    omp_event_handle_t event{};
    // where it started, in a trace
    double traceBegin{-1};
  };

  mutable std::mutex mutex_;
  std::vector<Exchange> handedOver_;
  // only the thread that polls touches these
  std::vector<Exchange> inFlight_;
  std::atomic<std::size_t> inFlightCount_{0};
  std::atomic<std::size_t> outstanding_{0};
};

/**
 * Runs the exchanges of a transport in the task graph of the host. A send becomes a detached task
 * that waits for the work that writes the copy layer before it, and that holds back the work that
 * overwrites the copy layer until the send has completed; likewise a receive, which waits for the
 * work that reads the ghost layer before it, and holds back the work that reads the ghost layer
 * after it until the data is there. So the exchanges count as launched as soon as they are started,
 * as on a device stream; the engine starts and tests them once their tasks are due.
 */
class TaskHaloTransport : public HaloTransport {
  public:
  TaskHaloTransport(std::unique_ptr<HaloTransport> transport, TaskHaloEngine& engine);
  ~TaskHaloTransport() override = default;

  TaskHaloTransport(const TaskHaloTransport&) = delete;
  TaskHaloTransport(TaskHaloTransport&&) = delete;
  TaskHaloTransport& operator=(const TaskHaloTransport&) = delete;
  TaskHaloTransport& operator=(TaskHaloTransport&&) = delete;

  void startSend() override;
  bool testSend() override { return true; }
  void startReceive() override;
  bool testReceive() override { return true; }
  [[nodiscard]] bool streamOrdered() const override { return true; }
  [[nodiscard]] std::string describe() const override;
  void startInterval(const ExchangeInterval& interval) override {
    transport_->startInterval(interval);
  }
  void finalize() override { transport_->finalize(); }
  void setHostTokens(const char* copy, char* ghost) override {
    copy_ = copy;
    ghost_ = ghost;
  }
  HaloTransport& dataTransport() override { return transport_->dataTransport(); }
  void setTraceName(const std::string& name) override { traceName_ = name; }

  /// records an exchange that began at `begin` (negative if the trace did not record then)
  void traceExchange(bool send, double begin) const;

  private:
  std::unique_ptr<HaloTransport> transport_;
  TaskHaloEngine& engine_;
  const char* copy_{nullptr};
  char* ghost_{nullptr};
  std::string traceName_{"halo exchange"};
};

} // namespace seissol::solver

#else

namespace seissol::solver {
/// without OpenMP 5.0, there is no task graph to run the exchanges in
class TaskHaloEngine {};
} // namespace seissol::solver

#endif // SEISSOL_HOST_TASK_GRAPH

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_TASKHALOTRANSPORT_H_
