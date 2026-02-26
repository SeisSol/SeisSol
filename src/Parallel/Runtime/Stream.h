// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_

#include "Memory/Tree/Layer.h"

#include <functional>
#include <utility>

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

  ManagedStream(ManagedStream&& old) noexcept : streamPtr_(old.streamPtr_) {}

  auto operator=(ManagedStream&& old) noexcept -> ManagedStream& {
    this->streamPtr_ = old.streamPtr_;
    return *this;
  }

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

  ManagedEvent(ManagedEvent&& old) noexcept : eventPtr_(old.eventPtr_) {}

  auto operator=(ManagedEvent&& old) noexcept -> ManagedEvent& {
    this->eventPtr_ = old.eventPtr_;
    return *this;
  }

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

  StreamRuntime() : StreamRuntime(0) {}

  explicit StreamRuntime(size_t ringbufferSize)
      : ringbufferSize_(ringbufferSize), disposed_(false) {
    streamPtr_.emplace();
    ringbufferPtr_.resize(ringbufferSize);
    events_.resize(EventPoolSize);

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
    device().api->streamHostFunction(streamPtr_->get(), std::forward<F>(handler));
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
    const auto pos = this->eventpos_;
    this->eventpos_ = (this->eventpos_ + 1) % this->events_.size();
    return this->events_[pos].get();
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
