// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_EVENTPOOL_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_EVENTPOOL_H_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace seissol::parallel::runtime {

namespace internal {
struct EventSlot {
  void* event{nullptr};
  std::atomic<std::uint32_t> refs{0};
  // the number of referenced slots of the pool the slot belongs to
  std::atomic<std::size_t>* referenced{nullptr};
};
} // namespace internal

/**
 * Reference to an event from an EventPool.
 *
 * The pool hands out a slot only while nobody refers to it, so an event cannot be re-recorded
 * underneath code that still intends to wait on it. Completion is not the criterion and cannot
 * be: an event may well have been reached and still be needed, and asking the device whether it
 * has been reached is not even allowed while a graph is being recorded - the query invalidates
 * the capture.
 *
 * Hold a reference from recording the event until the wait on it has been enqueued. After that
 * the wait no longer depends on the event, because both CUDA and HIP take the event's state at
 * the time the wait is issued.
 *
 * References may be copied and dropped on any thread.
 */
class EventRef {
  public:
  EventRef() = default;
  explicit EventRef(internal::EventSlot* slot) : slot_(slot) { acquire(); }

  EventRef(const EventRef& other) : slot_(other.slot_) { acquire(); }
  EventRef(EventRef&& other) noexcept : slot_(std::exchange(other.slot_, nullptr)) {}

  auto operator=(const EventRef& other) -> EventRef& {
    if (this != &other) {
      release();
      slot_ = other.slot_;
      acquire();
    }
    return *this;
  }

  auto operator=(EventRef&& other) noexcept -> EventRef& {
    if (this != &other) {
      release();
      slot_ = std::exchange(other.slot_, nullptr);
    }
    return *this;
  }

  ~EventRef() { release(); }

  [[nodiscard]] bool isValid() const { return slot_ != nullptr; }

  [[nodiscard]] void* get() const { return slot_ == nullptr ? nullptr : slot_->event; }

  private:
  void acquire() {
    if (slot_ != nullptr && slot_->refs.fetch_add(1, std::memory_order_relaxed) == 0) {
      slot_->referenced->fetch_add(1, std::memory_order_relaxed);
    }
  }

  void release() {
    if (slot_ != nullptr) {
      if (slot_->refs.fetch_sub(1, std::memory_order_release) == 1) {
        slot_->referenced->fetch_sub(1, std::memory_order_relaxed);
      }
      slot_ = nullptr;
    }
  }

  internal::EventSlot* slot_{nullptr};
};

/**
 * A pool of events that recycles an event only once nobody refers to it any more (see EventRef).
 * It grows when all of its events are referenced.
 *
 * The events come from the device by default. For tests, any other way to create and destroy them
 * can be plugged in.
 *
 * Handing out events is for one thread at a time; the references may live anywhere.
 */
class EventPool {
  public:
  using Create = std::function<void*()>;
  using Destroy = std::function<void(void*)>;

  static constexpr std::size_t InitialSize = 100;
  static constexpr std::size_t Growth = 100;
  static constexpr std::size_t MaxSize = 100000;

  /**
   * A pool of device events; without a device, of null events.
   */
  EventPool();

  EventPool(Create create, Destroy destroy, std::size_t initialSize = InitialSize);

  ~EventPool();

  EventPool(const EventPool&) = delete;
  EventPool(EventPool&&) = delete;
  EventPool& operator=(const EventPool&) = delete;
  EventPool& operator=(EventPool&&) = delete;

  /**
   * An event that nobody refers to.
   */
  EventRef next();

  /**
   * Destroys the events; needs to happen before the device is finalized. An event that is still
   * referenced is reported and left alone, neither destroyed nor freed: somebody may still wait
   * for it. Returns the number of such events.
   */
  std::size_t dispose();

  /// the number of events of the pool
  [[nodiscard]] std::size_t size() const { return slots_.size(); }

  /// the number of events that are referenced right now
  [[nodiscard]] std::size_t referenced() const {
    return referenced_ == nullptr ? 0 : referenced_->load(std::memory_order_relaxed);
  }

  /// the largest number of events that were referenced at once when an event was handed out
  [[nodiscard]] std::size_t peakReferenced() const { return peakReferenced_; }

  private:
  void grow(std::size_t count);

  Create create_;
  Destroy destroy_;
  // slots are held indirectly so that their addresses survive the pool growing; an EventRef
  // points straight at its slot, and so does it at the count of referenced slots
  std::vector<std::unique_ptr<internal::EventSlot>> slots_;
  std::unique_ptr<std::atomic<std::size_t>> referenced_;
  std::size_t position_{0};
  std::size_t peakReferenced_{0};
  bool disposed_{false};
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_EVENTPOOL_H_
