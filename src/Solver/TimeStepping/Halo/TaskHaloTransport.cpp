// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "TaskHaloTransport.h"

#include "Parallel/Runtime/HostRuntime.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <atomic>
#include <cstddef>
#include <memory>
#include <mutex>
#include <omp.h>
#include <sstream>
#include <string>
#include <utility>

#ifdef SEISSOL_HOST_TASK_GRAPH

namespace seissol::solver {

void TaskHaloEngine::enqueue(HaloTransport* transport, bool send, omp_event_handle_t event) {
  const std::scoped_lock lock(mutex_);
  handedOver_.push_back({transport, send, event});
}

void TaskHaloEngine::poll() {
  {
    const std::scoped_lock lock(mutex_);
    for (const auto& exchange : handedOver_) {
      if (exchange.send) {
        exchange.transport->startSend();
      } else {
        exchange.transport->startReceive();
      }
      inFlight_.push_back(exchange);
    }
    handedOver_.clear();
  }
  for (std::size_t i = 0; i < inFlight_.size();) {
    const auto& exchange = inFlight_[i];
    const bool completed =
        exchange.send ? exchange.transport->testSend() : exchange.transport->testReceive();
    if (completed) {
      omp_fulfill_event(exchange.event);
      outstanding_.fetch_sub(1, std::memory_order_release);
      inFlight_[i] = inFlight_.back();
      inFlight_.pop_back();
    } else {
      ++i;
    }
  }
  inFlightCount_.store(inFlight_.size(), std::memory_order_relaxed);
}

std::string TaskHaloEngine::describe() const {
  std::size_t handedOver = 0;
  {
    const std::scoped_lock lock(mutex_);
    handedOver = handedOver_.size();
  }
  std::ostringstream text;
  text << outstanding() << " tasks of exchanges in the task graph, " << handedOver
       << " of them ready to start, " << inFlightCount_.load(std::memory_order_relaxed)
       << " in flight";
  return text.str();
}

TaskHaloTransport::TaskHaloTransport(std::unique_ptr<HaloTransport> transport,
                                     TaskHaloEngine& engine)
    : transport_(std::move(transport)), engine_(engine) {}

void TaskHaloTransport::startSend() {
  HaloTransport* transport = transport_.get();
  TaskHaloEngine* engine = &engine_;
  // the detach clause sets the event when it creates the task
  // NOLINTNEXTLINE(misc-const-correctness,bugprone-invalid-enum-default-initialization)
  omp_event_handle_t event{};
  engine_.created();
#pragma omp task firstprivate(transport, engine) detach(event) depend(in : this -> copy_[0])
  engine->enqueue(transport, true, event);
}

void TaskHaloTransport::startReceive() {
  HaloTransport* transport = transport_.get();
  TaskHaloEngine* engine = &engine_;
  // the detach clause sets the event when it creates the task
  // NOLINTNEXTLINE(misc-const-correctness,bugprone-invalid-enum-default-initialization)
  omp_event_handle_t event{};
  engine_.created();
#pragma omp task firstprivate(transport, engine) detach(event) depend(out : this -> ghost_[0])
  engine->enqueue(transport, false, event);
}

std::string TaskHaloTransport::describe() const {
  const auto inner = transport_->describe();
  return inner.empty() ? engine_.describe() : inner + "; " + engine_.describe();
}

} // namespace seissol::solver

#endif // SEISSOL_HOST_TASK_GRAPH
