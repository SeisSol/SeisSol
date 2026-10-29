// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_EVENT_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_EVENT_H_

#include "Parallel/Runtime/EventPool.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace seissol::parallel::runtime {

/**
 * A counter in host memory that the device can read (see device::AbstractAPI::streamWaitMemory()),
 * and the value to wait for.
 */
struct HostCounter {
  const std::uint32_t* counter{nullptr};
  std::uint32_t value{0};

  [[nodiscard]] bool operator==(const HostCounter& other) const {
    return counter == other.counter && value == other.value;
  }
};

/**
 * Completes with the work enqueued up to where it was taken (see StreamRuntime::eventRecord()), so
 * that other work can be ordered after it (see StreamRuntime::eventSync()).
 *
 * Running directly, it is an event recorded on a stream. One from an event pool stays reserved
 * while held: the pool does not hand it out to be recorded anew (see EventRef). One that does not
 * come from a pool is not reserved; whoever hands it out has to keep it valid, and must not record
 * it anew, while anybody may wait for it.
 *
 * While the work gets built into a graph node by node, it is the node the work ends with instead;
 * the work that waits for it gets an edge from that node. Only work of the same graph can wait for
 * a node.
 *
 * Or it stands for counters in host memory that the device can read, which the host advances (e.g.
 * once MPI has completed an operation, see parallel::runtime::HostActor): it completes once each of
 * them has reached its value.
 */
class Event {
  public:
  using HostCounter = runtime::HostCounter;

  Event() = default;

  /// an event that does not come from a pool; not reserved
  explicit Event(void* event) : event_(event) {}

  /// an event from an event pool; reserved while held
  explicit Event(EventRef event) : event_(event.get()), reference_(std::move(event)) {}

  /// a node of a graph that is being built: the graph (as the device names it), and the node in it
  Event(const void* graph, std::size_t node) : graph_(graph), node_(node) {}

  /// counters in host memory, each of which has to reach its value; the second one may be empty
  explicit Event(HostCounter first, HostCounter second = HostCounter{})
      : counters_{first, second} {}

  /// the device event; null for a node
  [[nodiscard]] void* get() const { return event_; }

  /// the graph of a node; null for a device event
  [[nodiscard]] const void* graph() const { return graph_; }

  /// the node in its graph
  [[nodiscard]] std::size_t node() const { return node_; }

  /// the counters in host memory; empty ones have no counter
  [[nodiscard]] const std::array<HostCounter, 2>& counters() const { return counters_; }

  /// whether it stands for counters in host memory
  [[nodiscard]] bool onCounters() const { return counters_[0].counter != nullptr; }

  /// whether there is anything to wait for
  explicit operator bool() const { return event_ != nullptr || graph_ != nullptr || onCounters(); }

  /// whether both stand for the same event, the same node, or the same counters
  [[nodiscard]] bool operator==(const Event& other) const {
    return event_ == other.event_ && graph_ == other.graph_ && node_ == other.node_ &&
           counters_ == other.counters_;
  }

  [[nodiscard]] bool operator!=(const Event& other) const { return !(*this == other); }

  private:
  void* event_{nullptr};
  EventRef reference_;
  const void* graph_{nullptr};
  std::size_t node_{0};
  std::array<HostCounter, 2> counters_{};
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_EVENT_H_
