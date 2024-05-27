#ifndef SEISSOL_PARALLEL_RUNTIME_STREAM_HPP
#define SEISSOL_PARALLEL_RUNTIME_STREAM_HPP

#include <Initializer/tree/Layer.hpp>
#include <functional>
#include <omp.h>
#include <utility>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::parallel::runtime {

enum class Runtime { Native, Sycl, OpenMP };

class StreamRuntime {
#ifdef ACL_DEVICE
  private:
  static device::DeviceInstance& device() { return device::DeviceInstance::getInstance(); }

  public:
  static constexpr size_t RingbufferSize = 4;

  StreamRuntime() : disposed(false) {
    streamPtr = device().api->createGenericStream();
    ringbufferPtr.resize(RingbufferSize);
    for (size_t i = 0; i < RingbufferSize; ++i) {
      ringbufferPtr[i] = device().api->createGenericStream();
    }

    allStreams.resize(RingbufferSize + 1);
    allStreams[0] = streamPtr;
    for (size_t i = 0; i < RingbufferSize; ++i) {
      allStreams[i + 1] = ringbufferPtr[i];
    }

    for (size_t i = 0; i < 2 * RingbufferSize + 2; ++i) {
      events.push_back(device().api->createEvent());
    }
  }

  void dispose() {
    if (!disposed) {
      device().api->destroyGenericStream(streamPtr);
      for (size_t i = 0; i < RingbufferSize; ++i) {
        device().api->destroyGenericStream(ringbufferPtr[i]);
      }
      for (auto* event : events) {
        device().api->destroyEvent(event);
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
  void runGraph(seissol::initializer::GraphKey computeGraphKey,
                seissol::initializer::Layer& layer,
                F&& handler) {
    auto computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

    if (!computeGraphHandle) {
      device().api->streamBeginCapture(allStreams);

      std::invoke(std::forward<F>(handler), *this);

      device().api->streamEndCapture();

      computeGraphHandle = device().api->getLastGraphHandle();
      layer.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    }

    if (computeGraphHandle.isInitialized()) {
      device().api->launchGraph(computeGraphHandle, streamPtr);
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

  void* getEvent() {
    for (auto* event : events) {
      if (device().api->isEventCompleted(event)) {
        return event;
      }
    }
    void* newEvent = device().api->createEvent();
    events.push_back(newEvent);
    return newEvent;
  }

  void dependency(StreamRuntime& other) {
    void* event = getEvent();
    device().api->recordEventOnStream(event, other.stream());
    device().api->syncStreamWithEvent(stream(), event);
  }

  private:
  bool disposed;
  void* streamPtr;
  std::vector<void*> ringbufferPtr;
  std::vector<void*> allStreams;
  std::vector<void*> events;
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
    enqueueHost([=]() {
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < elemCount; ++i) {
        std::invoke(handler, i);
      }
    });
  }
#endif
};

} // namespace seissol::parallel::runtime

#endif
