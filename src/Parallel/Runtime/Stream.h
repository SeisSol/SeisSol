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
#include <omp.h>
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

  StreamRuntime() : disposed(false) {
    streamPtr = device().api->createStream();
    ringbufferPtr.resize(RingbufferSize);
    forkEvents.resize(RingbufferSize);
    joinEvents.resize(RingbufferSize);
    for (size_t i = 0; i < RingbufferSize; ++i) {
      ringbufferPtr[i] = device().api->createStream();
      forkEvents[i] = device().api->createEvent();
      joinEvents[i] = device().api->createEvent();
    }

    allStreams.resize(RingbufferSize + 1);
    allStreams[0] = streamPtr;
    for (size_t i = 0; i < RingbufferSize; ++i) {
      allStreams[i + 1] = ringbufferPtr[i];
    }

    forkEventSycl = device().api->createEvent();
    joinEventSycl = device().api->createEvent();
  }

  void dispose() {
    if (!disposed) {
      device().api->destroyGenericStream(streamPtr);
      for (size_t i = 0; i < RingbufferSize; ++i) {
        device().api->destroyGenericStream(ringbufferPtr[i]);
        device().api->destroyEvent(forkEvents[i]);
        device().api->destroyEvent(joinEvents[i]);
      }
      device().api->destroyEvent(forkEventSycl);
      device().api->destroyEvent(joinEventSycl);
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
  void enqueueOmpFor(std::size_t elemCount, F&& handler) {
    enqueueHost([=]() {
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < elemCount; ++i) {
        std::invoke(handler, i);
      }
    });
  }

  template <typename F>
  void envMany(size_t count, F&& handler) {
    for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
      device().api->recordEventOnStream(forkEvents[i], streamPtr);
      device().api->syncStreamWithEvent(ringbufferPtr[i], forkEvents[i]);
    }
    for (size_t i = 0; i < count; ++i) {
      std::invoke(handler, ringbufferPtr[i % ringbufferPtr.size()], i);
    }
    for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
      device().api->recordEventOnStream(joinEvents[i], ringbufferPtr[i]);
      device().api->syncStreamWithEvent(streamPtr, joinEvents[i]);
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

  template <typename F>
  void envSycl(void* queue, F&& handler) {
    syncToSycl(queue);
    std::invoke(handler);
    syncFromSycl(queue);
  }

  void syncToSycl(void* queue);
  void syncFromSycl(void* queue);

  /*
  // disabled unless using a modern compiler
    template <typename F>
    void envOMP(omp_depend_t& depobj, F&& handler) {
      syncToOMP(depobj);
      std::invoke(handler);
      syncFromOMP(depobj);
    }

    void syncToOMP(omp_depend_t& depobj);
    void syncFromOMP(omp_depend_t& depobj);
  */

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
  std::vector<void*> forkEvents;
  std::vector<void*> joinEvents;
  void* forkEventSycl;
  void* joinEventSycl;
#else
  public:
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
