// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "TaskHaloTransport.h"

#include "Monitoring/Trace.h"
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

void TaskHaloEngine::enqueue(const TaskHaloTransport* owner,
                             HaloTransport* transport,
                             bool send,
                             omp_event_handle_t event) {
  const std::scoped_lock lock(mutex_);
  handedOver_.push_back({owner, transport, send, event});
}

void TaskHaloEngine::poll() {
  {
    const std::scoped_lock lock(mutex_);
    for (auto exchange : handedOver_) {
      exchange.traceBegin = monitoring::Trace::recording() ? monitoring::Trace::now() : -1;
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
      exchange.owner->traceExchange(exchange.send, exchange.traceBegin);
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
  text << outstanding() << " exchanges outstanding (" << handedOver << " ready to start, "
       << inFlightCount_.load(std::memory_order_relaxed) << " in flight)";
  return text.str();
}

TaskHaloTransport::TaskHaloTransport(std::unique_ptr<HaloTransport> transport,
                                     TaskHaloEngine& engine)
    : transport_(std::move(transport)), engine_(engine) {}

void TaskHaloTransport::traceExchange(bool send, double begin) const {
  if (begin < 0 || !monitoring::Trace::recording()) {
    return;
  }
  static const auto SendName = monitoring::Trace::name("send");
  static const auto ReceiveName = monitoring::Trace::name("receive");
  // the tasks of the exchanges of a direction follow each other
  const auto track =
      send ? monitoring::Trace::track(static_cast<const void*>(&copy_), traceName_ + " send")
           : monitoring::Trace::track(static_cast<const void*>(&ghost_), traceName_ + " receive");
  monitoring::Trace::span(send ? SendName : ReceiveName, track, begin, monitoring::Trace::now());
}

void TaskHaloTransport::startSend() {
  const TaskHaloTransport* owner = this;
  HaloTransport* transport = transport_.get();
  TaskHaloEngine* engine = &engine_;
  // the detach clause sets the event when it creates the task
  // NOLINTNEXTLINE(misc-const-correctness,bugprone-invalid-enum-default-initialization)
  omp_event_handle_t event{};
  engine_.created();
#pragma omp task firstprivate(owner, transport, engine) detach(event) depend(in : this -> copy_[0])
  engine->enqueue(owner, transport, true, event);
}

void TaskHaloTransport::startReceive() {
  const TaskHaloTransport* owner = this;
  HaloTransport* transport = transport_.get();
  TaskHaloEngine* engine = &engine_;
  // the detach clause sets the event when it creates the task
  // NOLINTNEXTLINE(misc-const-correctness,bugprone-invalid-enum-default-initialization)
  omp_event_handle_t event{};
  engine_.created();
#pragma omp task firstprivate(owner, transport, engine) detach(event)                              \
    depend(out : this -> ghost_[0])
  engine->enqueue(owner, transport, false, event);
}

std::string TaskHaloTransport::describe() const {
  const auto inner = transport_->describe();
  return inner.empty() ? engine_.describe() : inner + "; " + engine_.describe();
}

} // namespace seissol::solver

#endif // SEISSOL_HOST_TASK_GRAPH
