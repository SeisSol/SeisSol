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
#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/EventPool.h"
#include "Parallel/Runtime/HostRuntime.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
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
 * super-timestep, see StreamRuntime::beginGraph()). All work then has to go through the runtimes
 * that take part in it; and no event may be recorded for timing, or queried.
 */
inline bool& recordingOuterGraph() {
  static thread_local bool recording = false;
  return recording;
}

/**
 * Whether host functions go onto the streams on HIP as well. Otherwise, HIP waits for the stream
 * and runs a host function right away. With concurrent clusters, the host must not wait for the
 * device during an action, though: the stream may wait for work of other clusters, or of other
 * processes (with an exchange ordered on the device), that the host enqueues only later.
 */
inline bool& hostFunctionsOnStreams() {
  static bool onStreams = false;
  return onStreams;
}

/**
 * Enqueues work on a stream of its own (and on further ones for independent branches, see
 * envMany()), or builds it into a graph.
 *
 * Runtimes form a hierarchy for graphs that span several of them: one runtime begins such a graph
 * (see beginGraph()), and the others attach to it (see attach()). Each of them then adds its work
 * as a chain of nodes, in the order it gets enqueued, and the work of one runtime waits for that of
 * another one where it says so (see eventRecord() and eventSync()): as edges of the graph, without
 * any event. The graphs of single actions (see runGraph()) become part of such a graph.
 */
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
      deferred_ = Event();
      streamPtr_.reset();
      ringbufferPtr_.clear();
      events_.dispose();
      disposed_ = true;
    }
  }

  ~StreamRuntime() { dispose(); }

  StreamRuntime(const StreamRuntime&) = delete;
  StreamRuntime& operator=(const StreamRuntime&) = delete;
  StreamRuntime(StreamRuntime&&) = delete;
  StreamRuntime& operator=(StreamRuntime&&) = delete;

  /**
   * Contributes the work of one callback to the current scope.
   *
   * While a graph is being built node by node, this becomes a single graph node that depends on
   * whatever was recorded before it. Otherwise it is just a call on the stream.
   */
  template <typename F>
  void record(F&& handler) {
    if (building()) {
      closeOpenNode();
      lastNodes_ = {addNode(takeDependencies(), streamPtr_->get(), std::forward<F>(handler))};
    } else {
      std::invoke(std::forward<F>(handler), mainStream());
    }
  }

  template <typename F>
  void env(F&& handler) {
    record(std::forward<F>(handler));
  }

  template <typename F>
  void enqueueHost(F&& handler) {
    // A graph that is being recorded cannot be waited for, but it takes the host function over.
    if (Backend != DeviceBackend::Hip || hostFunctionsOnStreams() || building() || capturing_ ||
        recordingOuterGraph()) {
      device().api().streamHostFunction(stream(), std::forward<F>(handler));
    } else {
      // if the stream host function call isn't implemented or slow, we'll need to synchronize
      wait();
      std::invoke(handler);
    }
  }

  template <typename F>
  void enqueueLoop(std::size_t elemCount,
                   const F& handler,
                   std::uint32_t traceName = HostRuntime::loopTraceName()) {
    enqueueHost([=]() { HostRuntime::forEach(elemCount, handler, traceName); });
  }

  EventRef nextEvent() { return events_.next(); }

  /**
   * Runs the handler for the branches `i` with `hasWork(i)` out of `count`, as envMany() does for
   * all of them; the others neither run nor get a node in a graph.
   */
  template <typename F, typename P>
  void envMany(size_t count, const F& handler, const P& hasWork) {
    std::vector<size_t> branches;
    for (size_t i = 0; i < count; ++i) {
      if (std::invoke(hasWork, i)) {
        branches.push_back(i);
      }
    }
    envMany(branches.size(),
            [&](void* stream, size_t branch) { std::invoke(handler, stream, branches[branch]); });
  }

  template <typename F>
  void envMany(size_t count, const F& handler) {
    if (count == 0) {
      return;
    }
    // a runtime attached to a graph without branches (see attach()) adds them to its chain
    if (building() && (parent_ == nullptr || branches_)) {
      closeOpenNode();

      // fork and join are properties of the graph here, so neither costs an event
      const auto forkDependencies = takeDependencies();

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

      // the work after it depends on all branches, which is exactly a join
      lastNodes_ = std::move(branches);
      return;
    }

    void* stream = building() ? this->stream() : mainStream();
    if (!building() && Backend != DeviceBackend::Hip && ringbufferSize_ > 0) {
      const auto forkEvent = nextEvent();
      device().api().recordEventOnStream(forkEvent.get(), stream);
      for (size_t i = 0; i < std::min(count, ringbufferPtr_.size()); ++i) {
        device().api().syncStreamWithEvent(ringbufferPtr_[i].get(), forkEvent.get());
      }
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, ringbufferPtr_[i % ringbufferPtr_.size()].get(), i);
      }
      for (size_t i = 0; i < std::min(count, ringbufferPtr_.size()); ++i) {
        const auto joinEvent = nextEvent();
        device().api().recordEventOnStream(joinEvent.get(), ringbufferPtr_[i].get());
        device().api().syncStreamWithEvent(stream, joinEvent.get());
      }
    } else {
      for (size_t i = 0; i < count; ++i) {
        std::invoke(handler, stream, i);
      }
    }
  }

  void wait() {
    if (building() || capturing_ || recordingOuterGraph()) {
      logError() << "Cannot synchronize a stream while a compute graph is being recorded on it.";
    }
    device().api().syncStreamWithHost(mainStream());
  }

  /**
   * The stream that work has to be issued on.
   *
   * Inside an explicit node this is that node's recording stream. While a graph is being built
   * from nodes and no explicit node is open, asking for the stream opens an implicit one, so
   * that code which was never taught about graphs still ends up in the graph rather than
   * running straight away. The implicit node is closed again at the next env, envMany, record,
   * eventRecord or eventSync, and at the end of the graph.
   *
   * The device captures one node of a graph at a time. So, while several runtimes build a graph
   * together, the implicit node also gets closed once another one of them starts a node: work must
   * not be issued on a stream taken from this one before another runtime has been called.
   */
  void* stream() {
    if (recordingStream_ != nullptr) {
      return recordingStream_;
    }
    if (building()) {
      openImplicitNode();
      return streamPtr_->get();
    }
    return mainStream();
  }

  template <typename F>
  void runGraphGeneric(device::DeviceGraphHandle& computeGraphHandle,
                       F&& handler,
                       GraphMode mode = GraphMode::Nodes) {
    if (building()) {
      // the work becomes part of the graph being built (e.g. that of a whole super-timestep); a
      // graph of its own can be neither built nor captured within it
      std::invoke(std::forward<F>(handler), *this);
      return;
    }
    // the graph gets built on the stream, which may not wait for anything while it does
    void* stream = mainStream();
    if (!computeGraphHandle.isInitialized()) {
      if (mode == GraphMode::Nodes && seissol::useGraphNodes()) {
        computeGraphHandle = device().api().graphCreate();

        buildGraph_ = computeGraphHandle;
        lastNodes_.clear();

        std::invoke(std::forward<F>(handler), *this);

        closeOpenNode();
        buildGraph_.reset();
        lastNodes_.clear();
        dependencies_.clear();
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
      device().api().launchGraph(computeGraphHandle, stream);
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

  /**
   * Begins a graph that spans this runtime and the ones that get attached to it (see attach()),
   * built node by node. The work that this runtime enqueues goes into it as well.
   */
  void beginGraph() {
    if (building() || capturing_) {
      logError() << "A stream runtime cannot begin a graph while it builds one already.";
    }
    buildGraph_ = device().api().graphCreate();
    if (!buildGraph_.isInitialized()) {
      logError() << "The device cannot build graphs node by node.";
    }
    lastNodes_.clear();
    dependencies_.clear();
  }

  /**
   * Ends the graph begun with beginGraph(): detaches all runtimes from it, and returns it ready to
   * be launched (see launchGraph()), as often as needed.
   */
  device::DeviceGraphHandle endGraph() {
    while (!children_.empty()) {
      children_.back()->detach();
    }
    closeImplicitNode();
    auto graph = std::exchange(buildGraph_, device::DeviceGraphHandle());
    lastNodes_.clear();
    dependencies_.clear();
    device().api().graphInstantiate(graph);
    return graph;
  }

  /**
   * Launches a graph (see endGraph()) on the stream of this runtime.
   */
  void launchGraph(const device::DeviceGraphHandle& graph) {
    device().api().launchGraph(graph, mainStream());
  }

  /**
   * Makes the work that this runtime enqueues from now on part of the graph that `parent` builds
   * (see beginGraph()), until it gets detached. The graph has to be launched after what this
   * runtime has enqueued before, and after a wait it defers (see deferWait()); whoever launches it
   * sees to that. Without `branches`, the independent branches of envMany() run one after another,
   * in the chain of nodes of the runtime.
   */
  void attach(StreamRuntime& parent, bool branches = true) {
    if (building() || capturing_) {
      logError() << "A stream runtime cannot attach to a graph while it builds one.";
    }
    deferred_ = Event();
    parent_ = &parent;
    branches_ = branches;
    parent.children_.push_back(this);
    lastNodes_.clear();
    dependencies_.clear();
  }

  /**
   * Captures each node of this runtime into a graph of its own, which then becomes a single node of
   * the graph being built (see device::AbstractAPI::graphBeginIsolatedNode): for libraries that
   * fail when they get captured into the same graph several times, such as NCCL.
   */
  void isolateNodes() { isolated_ = true; }

  /**
   * Lets the work of this runtime run directly again (see attach()).
   */
  void detach() {
    if (parent_ != nullptr) {
      closeImplicitNode();
      auto& siblings = parent_->children_;
      siblings.erase(std::remove(siblings.begin(), siblings.end(), this), siblings.end());
      parent_ = nullptr;
      lastNodes_.clear();
      dependencies_.clear();
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

  /**
   * Makes the work enqueued from now on wait for the event: for a node of the graph that the work
   * goes into, the next node depends on it; otherwise, the stream waits for the event.
   */
  void eventSync(const Event& event) {
    if (!event) {
      return;
    }
    if (building()) {
      if (event.graph() != buildGraph().get()) {
        logError()
            << "Work that gets built into a graph can only wait for the nodes of that graph.";
      }
      closeImplicitNode();
      dependencies_.emplace_back(event.node());
      return;
    }
    if (event.graph() != nullptr) {
      logError() << "Work that runs directly cannot wait for a node of a graph.";
    }
    device().api().syncStreamWithEvent(stream(), event.get());
  }

  /**
   * An event that completes with the work enqueued so far (see Event): while the work goes into a
   * graph node by node, the node it ends with (none before the first node).
   */
  Event eventRecord() {
    if (building()) {
      closeImplicitNode();
      if (!dependencies_.empty() || lastNodes_.size() > 1) {
        // a node that stands for the work so far and what it has been ordered after since: an
        // empty recorder yields a handle that stands for its own dependencies
        closeOpenNode();
        lastNodes_ = {addNode(takeDependencies(), streamPtr_->get(), [](void*) {})};
      }
      return lastNodes_.empty() ? Event()
                                : Event(buildGraph().get(), lastNodes_.front().getNodeId());
    }
    auto event = nextEvent();
    device().api().recordEventOnStream(event.get(), stream());
    return Event(std::move(event));
  }

  /**
   * Makes the work enqueued from now on wait for the event, like eventSync(); but the wait only
   * gets enqueued with that work. Nothing gets enqueued while the runtime has nothing to do: e.g.
   * while its work gets replayed from a recording, on the stream of the recording, which needs
   * none of its waits. Replaces a wait deferred before, which the event has to come after; an
   * empty event lets go of it.
   */
  void deferWait(Event event) { deferred_ = std::move(event); }

  /**
   * Whether the wait for the event is still deferred (see deferWait()), i.e. the runtime has
   * enqueued nothing since.
   */
  [[nodiscard]] bool defersWaitFor(const Event& event) const {
    return static_cast<bool>(event) && deferred_ == event;
  }

  /**
   * Enqueues a deferred wait (see deferWait()) right away, and lets go of its event.
   */
  void enqueueDeferredWait() {
    if (deferred_) {
      mainStream();
    }
  }

  [[nodiscard]] const EventPool& eventPool() const { return events_; }

  private:
  /**
   * The stream of the runtime, for work that runs directly: with the wait that has been deferred to
   * such work (see deferWait()) enqueued first.
   */
  void* mainStream() {
    if (recordingOuterGraph() && !capturing_) {
      logError()
          << "A stream runtime got work while the calling thread records a graph without it.";
    }
    if (deferred_) {
      device().api().syncStreamWithEvent(streamPtr_->get(), deferred_.get());
      deferred_ = Event();
    }
    return streamPtr_->get();
  }

  /// the runtime that holds the graph that the work goes into: the one attached to, or this one
  StreamRuntime& root() { return parent_ != nullptr ? *parent_ : *this; }

  /// the graph that the work goes into, node by node; empty if none
  device::DeviceGraphHandle& buildGraph() { return root().buildGraph_; }

  [[nodiscard]] bool building() const {
    return (parent_ != nullptr ? parent_->buildGraph_ : buildGraph_).isInitialized();
  }

  /// what the next node of this runtime depends on: its last nodes, and what it waits for since
  std::vector<device::DeviceGraphNodeHandle> takeDependencies() {
    auto dependencies = std::exchange(dependencies_, {});
    dependencies.insert(dependencies.end(), lastNodes_.begin(), lastNodes_.end());
    return dependencies;
  }

  template <typename F>
  device::DeviceGraphNodeHandle
      addNode(const std::vector<device::DeviceGraphNodeHandle>& dependencies,
              void* streamPtr,
              F&& handler) {
    void* previousRecordingStream = recordingStream_;
    recordingStream_ = streamPtr;
    beginNode(dependencies, streamPtr);
    std::invoke(std::forward<F>(handler), streamPtr);
    auto node = device().api().graphEndNode(buildGraph(), streamPtr);
    recordingStream_ = previousRecordingStream;
    return node;
  }

  void beginNode(const std::vector<device::DeviceGraphNodeHandle>& dependencies, void* streamPtr) {
    if (isolated_) {
      device().api().graphBeginIsolatedNode(buildGraph(), dependencies, streamPtr);
    } else {
      device().api().graphBeginNode(buildGraph(), dependencies, streamPtr);
    }
  }

  /// closes the implicit node that is open in the graph, of whichever runtime builds it
  void closeOpenNode() {
    auto* open = root().openNode_;
    if (open != nullptr) {
      open->closeImplicitNode();
    }
  }

  void openImplicitNode() {
    auto& root = this->root();
    if (root.openNode_ != this) {
      closeOpenNode();
      beginNode(takeDependencies(), streamPtr_->get());
      root.openNode_ = this;
    }
  }

  void closeImplicitNode() {
    auto& root = this->root();
    if (root.openNode_ == this) {
      lastNodes_ = {device().api().graphEndNode(root.buildGraph_, streamPtr_->get())};
      root.openNode_ = nullptr;
    }
  }

  // the graph this runtime builds (see runGraphGeneric() and beginGraph()); and of that, the
  // runtime whose implicit node is open, and the runtimes attached to it
  device::DeviceGraphHandle buildGraph_;
  StreamRuntime* openNode_{nullptr};
  std::vector<StreamRuntime*> children_;
  // the runtime whose graph this one builds as well, and whether with branches (see attach())
  StreamRuntime* parent_{nullptr};
  bool branches_{true};
  // the nodes the work of this runtime ends with (several after independent branches, see
  // envMany()), and what it waits for since (see eventSync())
  std::vector<device::DeviceGraphNodeHandle> lastNodes_;
  std::vector<device::DeviceGraphNodeHandle> dependencies_;
  void* recordingStream_{nullptr};
  bool capturing_{false};
  bool isolated_{false};
  Event deferred_;

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
  void enqueueLoop(std::size_t elemCount,
                   const F& handler,
                   std::uint32_t traceName = HostRuntime::loopTraceName()) {
    HostRuntime::forEach(elemCount, handler, traceName);
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

  void eventSync(const Event& /*event*/) {}
  Event eventRecord() { return {}; }
  void deferWait(const Event& /*event*/) {}
  [[nodiscard]] bool defersWaitFor(const Event& /*event*/) const { return false; }
  void enqueueDeferredWait() {}
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
