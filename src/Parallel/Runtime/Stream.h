// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_

#include "Memory/Tree/Layer.h"

#include <algorithm>
#include <functional>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel::runtime {

template <typename T>
class StreamMemoryHandle;

class ManagedStream {
  public:
  ManagedStream();

  ManagedStream(const ManagedStream&) = delete;
  auto operator=(const ManagedStream&) = delete;

  ManagedStream(ManagedStream&& old) noexcept
      : streamPtr_(std::exchange(old.streamPtr_, nullptr)) {}

  auto operator=(ManagedStream&& old) noexcept -> ManagedStream&;

  [[nodiscard]] void* get() const { return streamPtr_; }

  ~ManagedStream();

  private:
  void* streamPtr_{nullptr};
};

class ManagedEvent {
  public:
  ManagedEvent();

  ManagedEvent(const ManagedEvent&) = delete;
  auto operator=(const ManagedEvent&) = delete;

  ManagedEvent(ManagedEvent&& old) noexcept : eventPtr_(std::exchange(old.eventPtr_, nullptr)) {}

  auto operator=(ManagedEvent&& old) noexcept -> ManagedEvent&;

  [[nodiscard]] void* get() const { return eventPtr_; }

  ~ManagedEvent();

  private:
  void* eventPtr_{nullptr};
};

class StreamRuntime {
#ifdef ACL_DEVICE
  private:
  static device::DeviceInstance& device() { return device::DeviceInstance::getInstance(); }

  public:
  static constexpr size_t EventPoolSize = 100;
  static constexpr size_t EventPoolGrowth = 100;
  static constexpr size_t MaxEventPoolSize = 100000;

  StreamRuntime() : StreamRuntime(0) {}

  explicit StreamRuntime(size_t ringbufferSize)
      : ringbufferSize_(ringbufferSize), disposed_(false) {
    streamPtr_.emplace();
    ringbufferPtr_.resize(ringbufferSize);
    events_.resize(EventPoolSize);
    recorded_.resize(EventPoolSize, false);

    allStreams_.resize(ringbufferSize + 1);
    allStreams_[0] = streamPtr_->get();
    for (size_t i = 0; i < ringbufferSize; ++i) {
      allStreams_[i + 1] = ringbufferPtr_[i].get();
    }
  }

  void dispose() {
    if (!disposed_) {
      streamPtr_.reset();
      ringbufferPtr_.clear();
      events_.clear();
      recorded_.clear();
      disposed_ = true;
    }
  }

  ~StreamRuntime() { dispose(); }

  template <typename F>
  void env(F&& handler) {
    std::invoke(std::forward<F>(handler), streamPtr_->get());
  }

  template <typename F>
  void enqueueHost(F&& handler) {
    if (Backend != DeviceBackend::Hip) {
      device().api->streamHostFunction(streamPtr_->get(), std::forward<F>(handler));
    } else {
      // if the stream host function call isn't implemented or slow, we'll need to synchronize
      wait();
      std::invoke(handler);
    }
  }

  template <typename F>
  void enqueueLoop(std::size_t elemCount, F&& handler) {
    enqueueHost([=]() {
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < elemCount; ++i) {
        std::invoke(handler, i);
      }
    });
  }

  void* nextEvent() {
    while (true) {
      const auto pos = this->eventpos_;
      void* event = this->events_[pos].get();

      // Handing out an event that has been recorded but not yet reached would make a later
      // wait observe a point in time that has already passed, silently and without any
      // ordering guarantee. Grow the pool instead of wrapping onto it.
      if (this->recorded_[pos] && !device().api->isEventCompleted(event)) {
        growEventPool();
        continue;
      }

      this->recorded_[pos] = true;
      this->eventpos_ = (this->eventpos_ + 1) % this->events_.size();
      return event;
    }
  }

  template <typename F>
  void envMany(size_t count, F&& handler) {
    if (Backend != DeviceBackend::Hip && ringbufferSize_ > 0) {
      void* forkEvent = nextEvent();
      device().api->recordEventOnStream(forkEvent, streamPtr_->get());
      for (size_t i = 0; i < std::min(count, ringbufferPtr_.size()); ++i) {
        device().api->syncStreamWithEvent(ringbufferPtr_[i].get(), forkEvent);
      }
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, ringbufferPtr_[i % ringbufferPtr_.size()].get(), i);
      }
      for (size_t i = 0; i < std::min(count, ringbufferPtr_.size()); ++i) {
        void* joinEvent = nextEvent();
        device().api->recordEventOnStream(joinEvent, ringbufferPtr_[i].get());
        device().api->syncStreamWithEvent(streamPtr_->get(), joinEvent);
      }
    } else {
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, streamPtr_->get(), i);
      }
    }
  }

  void wait() { device().api->syncStreamWithHost(streamPtr_->get()); }

  void* stream() { return streamPtr_->get(); }

  template <typename F>
  void runGraphGeneric(device::DeviceGraphHandle& computeGraphHandle, F&& handler) {
    if (!computeGraphHandle.isInitialized()) {
      computeGraphHandle = device().api->streamBeginCapture(allStreams_);

      std::invoke(std::forward<F>(handler), *this);

      device().api->streamEndCapture(computeGraphHandle);
    }

    if (computeGraphHandle.isInitialized()) {
      device().api->launchGraph(computeGraphHandle, streamPtr_->get());
    }
  }

  template <typename VarmapT, typename F>
  void runGraph(seissol::initializer::GraphKey computeGraphKey,
                initializer::Layer<VarmapT>& layer,
                F&& handler) {
    auto computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

    const bool needsUpdate = !computeGraphHandle.isInitialized();

    runGraphGeneric(computeGraphHandle, std::forward<F>(handler));

    if (needsUpdate && computeGraphHandle.isInitialized()) {
      layer.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    }
  }

  template <typename T>
  T* allocMemory(std::size_t count) {
    return reinterpret_cast<T*>(device().api->allocMemAsync(count * sizeof(T), streamPtr_->get()));
  }

  template <typename T>
  void freeMemory(T* ptr) {
    device().api->freeMemAsync(ptr, streamPtr_->get());
  }

  template <typename T>
  StreamMemoryHandle<T> memoryHandle(std::size_t count) {
    return StreamMemoryHandle<T>(count, *this);
  }

  void eventSync(void* event) { device().api->syncStreamWithEvent(streamPtr_->get(), event); }

  private:
  void growEventPool() {
    const auto oldSize = this->events_.size();
    if (oldSize >= MaxEventPoolSize) {
      logError() << "Ran out of device events (" << oldSize
                 << " in flight). This points at events being recorded but never waited "
                    "upon.";
    }

    // Appending keeps the handles that have already been handed out valid: the event handle
    // is stored by value, so moving a ManagedEvent during reallocation does not change it.
    this->events_.resize(std::min(oldSize + EventPoolGrowth, MaxEventPoolSize));
    this->recorded_.resize(this->events_.size(), false);
    this->eventpos_ = oldSize;
  }

  public:
  void* eventRecord() {
    void* event = nextEvent();
    device().api->recordEventOnStream(event, streamPtr_->get());
    return event;
  }

  private:
  std::size_t ringbufferSize_{0};
  bool disposed_;
  std::optional<ManagedStream> streamPtr_;
  std::vector<ManagedStream> ringbufferPtr_;
  std::vector<void*> allStreams_;
  std::vector<ManagedEvent> events_;
  std::vector<bool> recorded_;
  std::size_t eventpos_{0};
#else
  public:
  StreamRuntime() : StreamRuntime(0) {}

  explicit StreamRuntime(std::size_t ringbufferSize) {}

  template <typename F>
  void enqueueLoop(std::size_t elemCount, const F& handler) {
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < elemCount; ++i) {
      std::invoke(handler, i);
    }
  }

  template <typename F>
  void enqueueHost(F&& handler) {
    std::invoke(std::forward<F>(handler));
  }

  void* stream() {
    // dummy
    return nullptr;
  }

  template <typename T>
  T* allocMemory(std::size_t count) {
    return new T[count];
  }
  template <typename T>
  void freeMemory(T* ptr) {
    delete[] ptr;
  }
  template <typename T>
  StreamMemoryHandle<T> memoryHandle(std::size_t count) {
    return StreamMemoryHandle<T>(count, *this);
  }

  void wait() {}
  void dispose() {}

  void eventSync(void* event) {}
  void* eventRecord() { return nullptr; }
#endif
};

template <typename T>
class StreamMemoryHandle {
  public:
  StreamMemoryHandle(std::size_t count, StreamRuntime& runtime)
      : data_(runtime.allocMemory<T>(count)), runtime_(runtime) {}

  StreamMemoryHandle(const StreamMemoryHandle&) = delete;
  auto operator=(const StreamMemoryHandle& stream) = delete;

  StreamMemoryHandle(StreamMemoryHandle&&) = default;
  auto operator=(StreamMemoryHandle&& stream) -> StreamMemoryHandle& = default;

  T* get() { return data_; }

  [[nodiscard]] const T* get() const { return data_; }

  ~StreamMemoryHandle() { runtime_.freeMemory(data_); }

  private:
  T* data_;
  StreamRuntime& runtime_;
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
