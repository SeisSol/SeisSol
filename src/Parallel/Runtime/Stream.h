// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_

#include <Memory/Tree/Layer.h>
#include <functional>
#include <utility>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::parallel::runtime {

enum class Runtime { Native, Sycl, OpenMP };

template <typename T>
class StreamMemoryHandle;

class StreamRuntime {
#ifdef ACL_DEVICE
  private:
  static device::DeviceInstance& device() { return device::DeviceInstance::getInstance(); }

  public:
  static constexpr size_t RingbufferSize = 4;
  static constexpr size_t EventPoolSize = 100;

  StreamRuntime() : disposed(false) {
    streamPtr = device().api->createStream();
    ringbufferPtr.resize(RingbufferSize);
    events.resize(EventPoolSize);
    for (size_t i = 0; i < RingbufferSize; ++i) {
      ringbufferPtr[i] = device().api->createStream();
    }
    for (size_t i = 0; i < events.size(); ++i) {
      events[i] = device().api->createEvent();
    }

    allStreams.resize(RingbufferSize + 1);
    allStreams[0] = streamPtr;
    for (size_t i = 0; i < RingbufferSize; ++i) {
      allStreams[i + 1] = ringbufferPtr[i];
    }
  }

  void dispose() {
    if (!disposed) {
      device().api->destroyGenericStream(streamPtr);
      for (size_t i = 0; i < RingbufferSize; ++i) {
        device().api->destroyGenericStream(ringbufferPtr[i]);
      }
      for (size_t i = 0; i < events.size(); ++i) {
        device().api->destroyEvent(events[i]);
      }
      disposed = true;
    }
  }

  ~StreamRuntime() { dispose(); }

  StreamRuntime(const StreamRuntime&) = delete;
  StreamRuntime operator=(const StreamRuntime&) = delete;

  template <typename F>
  void env(F&& handler) {
    std::invoke(std::forward<F>(handler), streamPtr);
  }

  template <typename F>
  void enqueueHost(F&& handler) {
    device().api->streamHostFunction(streamPtr, std::forward<F>(handler));
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
    return events[pos];
  }

  template <typename F>
  void envMany(size_t count, F&& handler) {
    if constexpr (Backend != DeviceBackend::Hip) {
      void* forkEvent = nextEvent();
      device().api->recordEventOnStream(forkEvent, streamPtr);
      for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
        device().api->syncStreamWithEvent(ringbufferPtr[i], forkEvent);
      }
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, ringbufferPtr[i % ringbufferPtr.size()], i);
      }
      for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
        void* joinEvent = nextEvent();
        device().api->recordEventOnStream(joinEvent, ringbufferPtr[i]);
        device().api->syncStreamWithEvent(streamPtr, joinEvent);
      }
    } else {
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, streamPtr, i);
      }
    }
  }

  void wait() { device().api->syncStreamWithHost(streamPtr); }

  void* stream() { return streamPtr; }

  template <typename F>
  void runGraphGeneric(device::DeviceGraphHandle& computeGraphHandle, F&& handler) {
    if (!computeGraphHandle.isInitialized()) {
      computeGraphHandle = device().api->streamBeginCapture(allStreams);

      std::invoke(std::forward<F>(handler), *this);

      device().api->streamEndCapture(computeGraphHandle);
    }

    if (computeGraphHandle.isInitialized()) {
      device().api->launchGraph(computeGraphHandle, streamPtr);
    }
  }

  template <typename F>
  void runGraph(seissol::initializer::GraphKey computeGraphKey,
                seissol::initializer::Layer& layer,
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
    return reinterpret_cast<T*>(device().api->allocMemAsync(count * sizeof(T), streamPtr));
  }

  template <typename T>
  void freeMemory(T* ptr) {
    device().api->freeMemAsync(reinterpret_cast<T*>(ptr), streamPtr);
  }

  template <typename T>
  StreamMemoryHandle<T> memoryHandle(std::size_t count) {
    return StreamMemoryHandle<T>(count, *this);
  }

  private:
  bool disposed;
  void* streamPtr;
  std::vector<void*> ringbufferPtr;
  std::vector<void*> allStreams;
  std::vector<void*> events;
  std::size_t eventpos{0};
#else
  public:
  template <typename F>
  void enqueueLoop(std::size_t elemCount, F handler) {
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < elemCount; ++i) {
      std::invoke(handler, i);
    }
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
  void wait() {}
  void dispose() {}
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
