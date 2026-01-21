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

  ManagedStream(ManagedStream&& old) noexcept : streamPtr(old.streamPtr) {}

  auto operator=(ManagedStream&& old) noexcept -> ManagedStream& {
    this->streamPtr = old.streamPtr;
    return *this;
  }

  [[nodiscard]] void* get() const { return streamPtr; }

  ~ManagedStream();

  private:
  void* streamPtr{nullptr};
};

class ManagedEvent {
  public:
  ManagedEvent();

  ManagedEvent(const ManagedEvent&) = delete;
  auto operator=(const ManagedEvent&) = delete;

  ManagedEvent(ManagedEvent&& old) noexcept : eventPtr(old.eventPtr) {}

  auto operator=(ManagedEvent&& old) noexcept -> ManagedEvent& {
    this->eventPtr = old.eventPtr;
    return *this;
  }

  [[nodiscard]] void* get() const { return eventPtr; }

  ~ManagedEvent();

  private:
  void* eventPtr{nullptr};
};

class StreamRuntime {
#ifdef ACL_DEVICE
  private:
  static device::DeviceInstance& device() { return device::DeviceInstance::getInstance(); }

  public:
  static constexpr size_t EventPoolSize = 100;

  StreamRuntime() : StreamRuntime(0) {}

  explicit StreamRuntime(size_t ringbufferSize) : ringbufferSize(ringbufferSize), disposed(false) {
    streamPtr.emplace();
    ringbufferPtr.resize(ringbufferSize);
    events.resize(EventPoolSize);

    allStreams.resize(ringbufferSize + 1);
    allStreams[0] = streamPtr->get();
    for (size_t i = 0; i < ringbufferSize; ++i) {
      allStreams[i + 1] = ringbufferPtr[i].get();
    }
  }

  void dispose() {
    if (!disposed) {
      streamPtr.reset();
      ringbufferPtr.clear();
      events.clear();
      disposed = true;
    }
  }

  ~StreamRuntime() { dispose(); }

  template <typename F>
  void env(F&& handler) {
    std::invoke(std::forward<F>(handler), streamPtr->get());
  }

  template <typename F>
  void enqueueHost(F&& handler) {
    device().api->streamHostFunction(streamPtr->get(), std::forward<F>(handler));
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
    const auto pos = eventpos;
    eventpos = (eventpos + 1) % events.size();
    return events[pos].get();
  }

  template <typename F>
  void envMany(size_t count, F&& handler) {
    if (Backend != DeviceBackend::Hip && ringbufferSize > 0) {
      void* forkEvent = nextEvent();
      device().api->recordEventOnStream(forkEvent, streamPtr->get());
      for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
        device().api->syncStreamWithEvent(ringbufferPtr[i].get(), forkEvent);
      }
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, ringbufferPtr[i % ringbufferPtr.size()].get(), i);
      }
      for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
        void* joinEvent = nextEvent();
        device().api->recordEventOnStream(joinEvent, ringbufferPtr[i].get());
        device().api->syncStreamWithEvent(streamPtr->get(), joinEvent);
      }
    } else {
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, streamPtr->get(), i);
      }
    }
  }

  void wait() { device().api->syncStreamWithHost(streamPtr->get()); }

  bool test() { return device().api->isStreamWorkDone(streamPtr->get()); }

  void* stream() { return streamPtr->get(); }

  template <typename F>
  void runGraphGeneric(device::DeviceGraphHandle& computeGraphHandle, F&& handler) {
    if (!computeGraphHandle.isInitialized()) {
      computeGraphHandle = device().api->streamBeginCapture(allStreams);

      std::invoke(std::forward<F>(handler), *this);

      device().api->streamEndCapture(computeGraphHandle);
    }

    if (computeGraphHandle.isInitialized()) {
      device().api->launchGraph(computeGraphHandle, streamPtr->get());
    }
  }

  template <typename VarmapT, typename F>
  void runGraph(seissol::initializer::GraphKey computeGraphKey,
                initializer::Layer<VarmapT>& layer,
                F&& handler) {
    auto computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

    bool needsUpdate = !computeGraphHandle.isInitialized();

    runGraphGeneric(computeGraphHandle, std::forward<F>(handler));

    if (needsUpdate && computeGraphHandle.isInitialized()) {
      layer.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    }
  }

  template <typename T>
  T* allocMemory(std::size_t count) {
    return reinterpret_cast<T*>(device().api->allocMemAsync(count * sizeof(T), streamPtr->get()));
  }

  template <typename T>
  void freeMemory(T* ptr) {
    device().api->freeMemAsync(reinterpret_cast<T*>(ptr), streamPtr->get());
  }

  template <typename T>
  StreamMemoryHandle<T> memoryHandle(std::size_t count) {
    return StreamMemoryHandle<T>(count, *this);
  }

  void eventSync(void* event) { device().api->syncStreamWithEvent(streamPtr->get(), event); }

  void* eventRecord() {
    void* event = nextEvent();
    device().api->recordEventOnStream(event, streamPtr->get());
    return event;
  }

  private:
  std::size_t ringbufferSize{0};
  bool disposed;
  std::optional<ManagedStream> streamPtr;
  std::vector<ManagedStream> ringbufferPtr;
  std::vector<void*> allStreams;
  std::vector<ManagedEvent> events;
  std::size_t eventpos{0};
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
  bool test() { return true; }
  void dispose() {}

  void eventSync(void* event) {}
  void* eventRecord() { return nullptr; }
#endif
};

template <typename T>
class StreamMemoryHandle {
  public:
  StreamMemoryHandle(std::size_t count, StreamRuntime& runtime)
      : data(runtime.allocMemory<T>(count)), runtime(runtime) {}

  StreamMemoryHandle(const StreamMemoryHandle&) = delete;
  auto operator=(const StreamMemoryHandle& stream) = delete;

  StreamMemoryHandle(StreamMemoryHandle&&) = default;
  auto operator=(StreamMemoryHandle&& stream) -> StreamMemoryHandle& = default;

  T* get() { return data; }

  const T* get() const { return data; }

  ~StreamMemoryHandle() { runtime.freeMemory(data); }

  private:
  T* data;
  StreamRuntime& runtime;
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
