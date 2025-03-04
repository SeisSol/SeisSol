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
#include <functional>
#include <mutex>
#include <omp.h>
#include <queue>
#include <utility>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::parallel::runtime {

class StreamRuntime {
  private:
  std::shared_ptr<seissol::parallel::host::CpuExecutor> cpu;
#ifdef ACL_DEVICE
  private:
  static device::DeviceInstance& device() { return device::DeviceInstance::getInstance(); }
  static std::mutex mutexCPU;

  public:
  static constexpr size_t RingbufferSize = 4;

  StreamRuntime(const std::shared_ptr<seissol::parallel::host::CpuExecutor>& cpu,
                double priority = 0.0)
      : cpu(cpu), disposed(false), priority(priority) {
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
    for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
      void* event = getEvent();
      device().api->recordEventOnStream(event, streamPtr);
      device().api->syncStreamWithEvent(ringbufferPtr[i], event);
    }
    for (size_t i = 0; i < count; ++i) {
      std::invoke(handler, ringbufferPtr[i % ringbufferPtr.size()], i);
    }
    for (size_t i = 0; i < std::min(count, ringbufferPtr.size()); ++i) {
      void* event = getEvent();
      device().api->recordEventOnStream(event, ringbufferPtr[i]);
      device().api->syncStreamWithEvent(streamPtr, event);
    }
  }

  void wait() { device().api->syncStreamWithHost(streamPtr); }

  void* stream() { return streamPtr; }

  template <typename F>
  void runGraphGeneric(device::DeviceGraphHandle& computeGraphHandle, F&& handler) {
    if (!computeGraphHandle.isInitialized()) {
      device().api->streamBeginCapture(allStreams);

      std::invoke(std::forward<F>(handler), *this);

      device().api->streamEndCapture();

      computeGraphHandle = device().api->getLastGraphHandle();
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

  private:
  bool disposed;
  double priority;
  void* streamPtr;
  std::vector<void*> ringbufferPtr;
  std::vector<void*> allStreams;
  std::queue<void*> events;
#else
  public:
  void wait() {}
  void dispose() {}

  template <typename F>
  void enqueueHost(F&& handler) {
    std::invoke(std::forward<F>(handler));
  }

  template <typename F>
  void enqueueOmpFor(std::size_t elemCount, F&& handler) {
    if (elemCount > 0) {
      enqueueHost([=]() {
#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < elemCount; ++i) {
          std::invoke(handler, i);
        }
      });
    }
  }
#endif
};

class HostRuntime {
  public:
};

class NativeRuntime {
  public:
};

class SyclRuntime {
  public:
};

class OpenMPRuntime {
  public:
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
