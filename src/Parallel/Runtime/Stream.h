#ifndef SEISSOL_PARALLEL_RUNTIME_STREAM_HPP
#define SEISSOL_PARALLEL_RUNTIME_STREAM_HPP

#include <Initializer/Tree/Layer.h>
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
  static constexpr size_t EventbufferSize = 100; // >= 16 needed at the moment

  StreamRuntime() : disposed(false) {
    streamPtr = device().api->createGenericStream();
    ringbufferPtr.resize(RingbufferSize);
    events.resize(EventbufferSize);
    for (size_t i = 0; i < RingbufferSize; ++i) {
      ringbufferPtr[i] = device().api->createGenericStream();
    }
    for (size_t i = 0; i < EventbufferSize; ++i) {
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
      for (size_t i = 0; i < EventbufferSize; ++i) {
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
#ifdef SEISSOL_KERNELS_HIP
    if (captureState) {
      const std::size_t base = captureStreams.size();
      for (size_t i = 0; i < std::min(count, RingbufferSize); ++i) {
        captureStreams.emplace_back(device().api->createGenericStream());
        void* event = getEvent();
        device().api->recordEventOnStream(event, streamPtr);
        device().api->syncStreamWithEvent(captureStreams[i + base], event);
      }
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, captureStreams[base + (i % RingbufferSize)], i);
      }
      for (size_t i = 0; i < std::min(count, RingbufferSize); ++i) {
        void* event = getEvent();
        device().api->recordEventOnStream(event, captureStreams[i + base]);
        device().api->syncStreamWithEvent(streamPtr, event);
      }
      return;
    }
#endif
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
      // TODO: remove captureState once ROCm is stable enough
      device().api->streamBeginCapture(allStreams);
      captureState = true;

      std::invoke(std::forward<F>(handler), *this);

      captureState = false;
      device().api->streamEndCapture();

      for (auto* stream : captureStreams) {
        device().api->destroyGenericStream(stream);
      }
      captureStreams.resize(0);

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

  void* getEvent() {
    void* event = events.at(eventPos);
    eventPos = (eventPos + 1) % events.size();
    return event;
  }

  private:
  bool disposed;
  void* streamPtr;
  std::vector<void*> ringbufferPtr;
  std::vector<void*> allStreams;
  std::vector<void*> events;
  std::size_t eventPos{0};
  bool captureState{false};
  std::vector<void*> captureStreams;
#else
  public:
  void wait() {}
  void dispose() {}
#endif
};

} // namespace seissol::parallel::runtime

#endif
