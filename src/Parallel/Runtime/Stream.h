// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_STREAM_H_

#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"
#include "Parallel/Runtime/EventPool.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel::runtime {

/**
 * How a compute graph is put together.
 *
 * Capture records a whole stream and lets the backend derive the dependency structure from the
 * events on it. Nodes states the structure directly, one node per env/envMany scope, so no
 * events are involved at all.
 *
 * Nodes does not require the call site to be routed: work that simply takes stream() opens an
 * implicit node and lands in the graph as well. env, envMany and record exist to give that
 * work a shape - a named node that other nodes can depend on, or a set of branches that run
 * independently.
 */
enum class GraphMode { Capture, Nodes };

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

/**
 * Whether the calling thread records a graph that spans several stream runtimes (e.g. a whole
 * super-timestep). The graphs of single actions then run their work directly, so that it becomes
 * part of it.
 */
inline bool& recordingOuterGraph() {
  static thread_local bool recording = false;
  return recording;
}

class StreamRuntime {
#ifdef ACL_DEVICE
  private:
  static device::DeviceInstance& device() { return device::DeviceInstance::instance(); }

  public:
  StreamRuntime() : StreamRuntime(0) {}

  explicit StreamRuntime(size_t ringbufferSize)
      : ringbufferSize_(ringbufferSize), disposed_(false) {
    streamPtr_.emplace();
    ringbufferPtr_.resize(ringbufferSize);

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
      events_.dispose();
      disposed_ = true;
    }
  }

  ~StreamRuntime() { dispose(); }

  /**
   * Contributes the work of one callback to the current scope.
   *
   * While a graph is being built node by node, this becomes a single graph node that depends on
   * whatever was recorded before it. Otherwise it is just a call on the stream.
   */
  template <typename F>
  void record(F&& handler) {
    if (buildGraph_.isInitialized()) {
      closeImplicitNode();

      std::vector<device::DeviceGraphNodeHandle> dependencies;
      if (lastNode_.isInitialized()) {
        dependencies.push_back(lastNode_);
      }
      lastNode_ = addNode(dependencies, streamPtr_->get(), std::forward<F>(handler));
    } else {
      std::invoke(std::forward<F>(handler), streamPtr_->get());
    }
  }

  template <typename F>
  void env(F&& handler) {
    record(std::forward<F>(handler));
  }

  template <typename F>
  void enqueueHost(F&& handler) {
    // A graph that is being recorded cannot be waited for, but it takes the host function over.
    // That includes a graph spanning several stream runtimes that the calling thread records.
    if (Backend != DeviceBackend::Hip || buildGraph_.isInitialized() || capturing_ ||
        recordingOuterGraph()) {
      device().api().streamHostFunction(stream(), std::forward<F>(handler));
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

  EventRef nextEvent() { return events_.next(); }

  template <typename F>
  void envMany(size_t count, F&& handler) {
    if (buildGraph_.isInitialized()) {
      closeImplicitNode();

      // fork and join are properties of the graph here, so neither costs an event
      std::vector<device::DeviceGraphNodeHandle> forkDependencies;
      if (lastNode_.isInitialized()) {
        forkDependencies.push_back(lastNode_);
      }

      // Sibling branches go onto different streams. The graph edges are what orders them, so
      // for CUDA and HIP the choice of stream makes no difference - but SYCL derives edges from
      // the recorded queue, and its queues are in order, so two branches sharing a queue would
      // be serialized. Branches beyond the ringbuffer size share as they do on the stream path.
      std::vector<device::DeviceGraphNodeHandle> branches;
      branches.reserve(count);
      for (size_t i = 0; i < count; ++i) {
        void* branchStream = ringbufferPtr_.empty()
                                 ? streamPtr_->get()
                                 : ringbufferPtr_[i % ringbufferPtr_.size()].get();
        branches.push_back(addNode(forkDependencies, branchStream, [&](void* streamPtr) {
          std::invoke(handler, streamPtr, i);
        }));
      }

      // an empty recorder yields a handle that stands for its own dependencies, which is
      // exactly a join
      lastNode_ = addNode(branches, streamPtr_->get(), [](void*) {});
      return;
    }

    if (Backend != DeviceBackend::Hip && ringbufferSize_ > 0) {
      const auto forkEvent = nextEvent();
      device().api().recordEventOnStream(forkEvent.get(), streamPtr_->get());
      for (size_t i = 0; i < std::min(count, ringbufferPtr_.size()); ++i) {
        device().api().syncStreamWithEvent(ringbufferPtr_[i].get(), forkEvent.get());
      }
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, ringbufferPtr_[i % ringbufferPtr_.size()].get(), i);
      }
      for (size_t i = 0; i < std::min(count, ringbufferPtr_.size()); ++i) {
        const auto joinEvent = nextEvent();
        device().api().recordEventOnStream(joinEvent.get(), ringbufferPtr_[i].get());
        device().api().syncStreamWithEvent(streamPtr_->get(), joinEvent.get());
      }
    } else {
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, streamPtr_->get(), i);
      }
    }
  }

  void wait() {
    if (buildGraph_.isInitialized() || capturing_ || recordingOuterGraph()) {
      logError() << "Cannot synchronize a stream while a compute graph is being recorded on it.";
    }
    device().api().syncStreamWithHost(streamPtr_->get());
  }

  /**
   * The stream that work has to be issued on.
   *
   * Inside an explicit node this is that node's recording stream. While a graph is being built
   * from nodes and no explicit node is open, asking for the stream opens an implicit one, so
   * that code which was never taught about graphs still ends up in the graph rather than
   * running straight away. The implicit node is closed again at the next env, envMany or record,
   * and at the end of the graph.
   */
  void* stream() {
    if (recordingStream_ != nullptr) {
      return recordingStream_;
    }
    if (buildGraph_.isInitialized()) {
      openImplicitNode();
    }
    return streamPtr_->get();
  }

  template <typename F>
  void runGraphGeneric(device::DeviceGraphHandle& computeGraphHandle,
                       F&& handler,
                       GraphMode mode = GraphMode::Nodes) {
    if (recordingOuterGraph()) {
      // the work becomes part of the graph being recorded; neither building nor capturing a graph
      // of its own may start on a stream that is being captured already
      std::invoke(std::forward<F>(handler), *this);
      return;
    }
    if (!computeGraphHandle.isInitialized()) {
      if (mode == GraphMode::Nodes && seissol::useGraphNodes()) {
        computeGraphHandle = device().api().graphCreate();

        buildGraph_ = computeGraphHandle;
        lastNode_ = device::DeviceGraphNodeHandle();

        std::invoke(std::forward<F>(handler), *this);

        closeImplicitNode();
        buildGraph_.reset();
        device().api().graphInstantiate(computeGraphHandle);
      } else {
        computeGraphHandle = device().api().streamBeginCapture(allStreams_);

        // a backend without graph support hands back an empty handle and captures nothing; the
        // handler then simply runs, and may synchronize like it does outside of graphs
        capturing_ = computeGraphHandle.isInitialized();
        std::invoke(std::forward<F>(handler), *this);
        capturing_ = false;

        if (computeGraphHandle.isInitialized()) {
          device().api().streamEndCapture(computeGraphHandle);
        }
      }
    }

    if (computeGraphHandle.isInitialized()) {
      device().api().launchGraph(computeGraphHandle, streamPtr_->get());
    }
  }

  /**
   * Runs the handler, replaying a cached graph for it where possible.
   *
   * `cacheable` states whether this invocation will recur. A graph key that only ever occurs
   * once - most importantly the truncated time step that every cluster takes right before a
   * synchronization point - would otherwise capture and instantiate a graph that is replayed
   * a single time and then kept forever, so the number of live graphs would grow with the
   * number of synchronization points.
   */
  template <typename VarmapT, typename F>
  void runGraph(seissol::initializer::GraphKey computeGraphKey,
                initializer::Layer<VarmapT>& layer,
                F&& handler,
                bool cacheable = true,
                GraphMode mode = GraphMode::Nodes) {
    if (!cacheable) {
      std::invoke(std::forward<F>(handler), *this);
      return;
    }

    auto computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

    const bool needsUpdate = !computeGraphHandle.isInitialized();

    runGraphGeneric(computeGraphHandle, std::forward<F>(handler), mode);

    if (needsUpdate && computeGraphHandle.isInitialized()) {
      layer.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    }
  }

  template <typename T>
  T* allocMemory(std::size_t count) {
    return reinterpret_cast<T*>(device().api().allocMemAsync(count * sizeof(T), stream()));
  }

  template <typename T>
  void freeMemory(T* ptr) {
    device().api().freeMemAsync(ptr, stream());
  }

  template <typename T>
  StreamMemoryHandle<T> memoryHandle(std::size_t count) {
    return StreamMemoryHandle<T>(count, *this);
  }

  void eventSync(const EventRef& event) {
    device().api().syncStreamWithEvent(stream(), event.get());
  }

  EventRef eventRecord() {
    auto event = nextEvent();
    device().api().recordEventOnStream(event.get(), stream());
    return event;
  }

  private:
  template <typename F>
  device::DeviceGraphNodeHandle
      addNode(const std::vector<device::DeviceGraphNodeHandle>& dependencies,
              void* streamPtr,
              F&& handler) {
    void* previousRecordingStream = recordingStream_;
    recordingStream_ = streamPtr;
    auto node =
        device().api().graphAddNode(buildGraph_, dependencies, streamPtr, [&](void* recorded) {
          std::invoke(std::forward<F>(handler), recorded);
        });
    recordingStream_ = previousRecordingStream;
    return node;
  }

  void openImplicitNode() {
    if (!implicitNodeOpen_) {
      std::vector<device::DeviceGraphNodeHandle> dependencies;
      if (lastNode_.isInitialized()) {
        dependencies.push_back(lastNode_);
      }
      device().api().graphBeginNode(buildGraph_, dependencies, streamPtr_->get());
      implicitNodeOpen_ = true;
    }
  }

  void closeImplicitNode() {
    if (implicitNodeOpen_) {
      lastNode_ = device().api().graphEndNode(buildGraph_, streamPtr_->get());
      implicitNodeOpen_ = false;
    }
  }

  device::DeviceGraphHandle buildGraph_;
  device::DeviceGraphNodeHandle lastNode_;
  void* recordingStream_{nullptr};
  bool implicitNodeOpen_{false};
  bool capturing_{false};

  std::size_t ringbufferSize_{0};
  bool disposed_;
  std::optional<ManagedStream> streamPtr_;
  std::vector<ManagedStream> ringbufferPtr_;
  std::vector<void*> allStreams_;
  EventPool events_;
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

  void eventSync(const EventRef& event) {}
  EventRef eventRecord() { return EventRef(); }
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
