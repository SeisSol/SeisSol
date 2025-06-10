// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_

#include <Memory/Tree/Layer.h>
#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Host/SyncExecutor.h>
#include <functional>
#include <mutex>
#include <omp.h>
#include <queue>
#include <utility>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::parallel::runtime {

#ifdef ACL_DEVICE
using EventT = void*;
#else
using EventT = std::vector<std::shared_ptr<host::Task>>;
#endif

template <typename T>
class StreamMemoryHandle;

class StreamRuntime {
  private:
  std::shared_ptr<seissol::parallel::host::CpuExecutor> cpu;
  static std::mutex mutexCPU;
#ifdef ACL_DEVICE
  private:
  static device::DeviceInstance& device() { return device::DeviceInstance::getInstance(); }

  public:
  static constexpr size_t RingbufferSize = 4;
  static constexpr size_t EventPoolSize = 100;

  StreamRuntime(const std::shared_ptr<seissol::parallel::host::CpuExecutor>& cpu =
                    std::make_shared<host::SyncExecutor>(),
                double priority = 0.0)
      : cpu(cpu), priority(priority) {
    streamPtr = device().api->createStream(priority);
    ringbufferPtr.resize(RingbufferSize);
    for (size_t i = 0; i < RingbufferSize; ++i) {
      ringbufferPtr[i] = device().api->createStream(priority);
    }

    allStreams.resize(RingbufferSize + 1);
    allStreams[0] = streamPtr;
    for (size_t i = 0; i < RingbufferSize; ++i) {
      allStreams[i + 1] = ringbufferPtr[i];
    }

    // heuristic value
    constexpr std::size_t StartEvents = 10000;
    for (size_t i = 0; i < StartEvents; ++i) {
      events.push(device().api->createEvent());
    }
  }

  void dispose() {
    if (!disposed) {
      device().api->destroyGenericStream(streamPtr);
      for (size_t i = 0; i < RingbufferSize; ++i) {
        device().api->destroyGenericStream(ringbufferPtr[i]);
      }
      while (!events.empty()) {
        device().api->destroyEvent(events.front());
        events.pop();
      }
    }
    disposed = true;
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
    if (elemCount > 0) {
      auto cpu = this->cpu;
      enqueueHost([cpu, elemCount, handler]() {
#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < elemCount; ++i) {
          std::invoke(handler, i);
        }
        // cpu->add(0, elemCount, handler, {})->wait();
      });
    }
  }

  template <typename F>
  void envMany(size_t count, F&& handler) {
    if constexpr (Backend != DeviceBackend::Hip) {
      void* forkEvent = getEvent();
      device().api->recordEventOnStream(forkEvent, streamPtr);
      for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
        device().api->syncStreamWithEvent(ringbufferPtr[i], forkEvent);
      }
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, ringbufferPtr[i % ringbufferPtr.size()], i);
      }
      for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
        void* joinEvent = getEvent();
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

  void recycleEvent(void* eventPtr) { events.push(eventPtr); }

  void* getEvent() {
    constexpr std::size_t NextEvents = 100;
    if (events.size() < NextEvents) {
      for (size_t i = 0; i < NextEvents; ++i) {
        events.push(device().api->createEvent());
      }
    }
    void* newEvent = events.front();
    events.pop();
    return newEvent;
  }

  void* recordEvent() {
    void* newEvent = getEvent();
    device().api->recordEventOnStream(newEvent, stream());
    return newEvent;
  }

  void waitEvent(void* eventPtr) { device().api->syncStreamWithEvent(stream(), eventPtr); }

  template <typename T>
  T* allocMemory(std::size_t count) {
    return reinterpret_cast<T*>(device().api->allocMemAsync(count * sizeof(T), streamPtr));
  }

  template <typename T>
  void freeMemory(T* ptr) {
    device().api->freeMemAsync(ptr, streamPtr);
  }

  template <typename T>
  StreamMemoryHandle<T> memoryHandle(std::size_t count) {
    return StreamMemoryHandle<T>(count, *this);
  }

  private:
  bool disposed{};
  double priority;
  void* streamPtr;
  std::vector<void*> ringbufferPtr;
  std::vector<void*> allStreams;
  std::queue<void*> events;
#else
  public:
  StreamRuntime(const std::shared_ptr<seissol::parallel::host::CpuExecutor>& cpu =
                    std::make_shared<host::SyncExecutor>(),
                double priority = 0.0)
      : cpu(cpu) {}
  void wait() {
    for (auto& task : waitTasks) {
      task->wait();
    }
  }
  template <typename T>
  T* allocMemory(std::size_t count) {
    return new T[count];
  }
  template <typename T>
  void freeMemory(T* ptr) {
    delete[] ptr;
  }
  void dispose() {}

  template <typename F>
  void enqueueHost(F&& handler) {
    auto last = cpu->add(0, 1, [handler](std::size_t) { std::invoke(handler); }, waitTasks);
    waitTasks = std::vector<std::shared_ptr<host::Task>>();
    waitTasks.emplace_back(last);
  }

  template <typename F>
  void enqueueOmpFor(std::size_t elemCount, F&& handler) {
    if (elemCount > 0) {
      auto last = cpu->add(0, elemCount, handler, waitTasks);
      waitTasks = std::vector<std::shared_ptr<host::Task>>();
      waitTasks.emplace_back(last);
    }
  }

  void* getEvent() { return nullptr; }

  EventT recordEvent() { return waitTasks; }

  void waitEvent(EventT eventPtr) {
    for (const auto& event : eventPtr) {
      waitTasks.emplace_back(event);
    }
  }

  void recycleEvent(EventT eventPtr) {
    // nop
  }

  private:
  std::vector<std::shared_ptr<host::Task>> waitTasks;
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
