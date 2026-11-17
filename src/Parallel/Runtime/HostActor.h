// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_HOSTACTOR_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_HOSTACTOR_H_

#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/HostBufferPool.h"
#include "Parallel/Runtime/HostRuntime.h"

#include <atomic>
#include <cstddef>
#include <deque>
#include <functional>
#include <memory>
#include <vector>

namespace seissol::parallel::runtime {

class StreamRuntime;

/**
 * Runs work on the host that reads what the device has written (e.g. the receiver samples, from
 * the DOFs the device has gathered into host memory) once the device is done with it, instead of a
 * host function on a stream. The thread that enqueues the work of the device polls the events the
 * work waits for (see progress()), and hands over what is ready. So the device never waits for the
 * host, and no host function costs its launch (up to half a millisecond on newer ROCm) or runs on
 * a thread of the device runtime.
 *
 * The work goes into queues: the work of one queue runs in the order it was posted, and the queues
 * do not wait for each other. Ready work becomes an action of the host runtime of its queue: with
 * a task graph on the host, a task after the work of the queue before it, which the other threads
 * of the epoch take; otherwise, it runs right away, with the threads of the polling thread in its
 * loops. Or an executor takes it instead (see setExecutor(); e.g. the thread of the outputs), which
 * runs what it gets in that order.
 *
 * Before the host reads what the work of a queue writes (e.g. the samples at a synchronization
 * point), drain() waits until all of it is done. Whatever the device writes for the work has to
 * stay untouched until the work has run: in buffers of its own, which the queues hand out (see
 * Queue::take()). The work may run on another thread, so all that it changes besides has to be its
 * own, or safe to change from there.
 *
 * While a graph gets recorded (e.g. a super-timestep, see beginRecording()), the work goes into the
 * recording instead, which posts it after each launch of the graph (see Recording).
 *
 * Only for the thread that enqueues the work of the device; no event gets queried while a graph is
 * being recorded (see recordingOuterGraph()).
 */
class HostActor {
  public:
  using Completed = std::function<bool(const Event&)>;
  using Executor = std::function<void(std::function<void()>)>;

  /**
   * Where the device writes the data of one piece of work (see Queue::take()): `device` for the
   * device, and `host` where the work reads it; null while a graph gets recorded, whose launches
   * copy it.
   */
  struct Staging {
    void* device{nullptr};
    void* host{nullptr};
  };

  class Recording;

  /**
   * The work of one producer on the device (e.g. the receivers of a time cluster), in the order
   * it was posted.
   */
  class Queue {
public:
    /// a queue of the actor of the process, whose pieces of work read `bytes` each (see take())
    explicit Queue(std::size_t bytes = 0);

    explicit Queue(HostActor& actor, std::size_t bytes = 0);

    ~Queue();

    Queue(const Queue&) = delete;
    Queue& operator=(const Queue&) = delete;
    Queue(Queue&&) = delete;
    Queue& operator=(Queue&&) = delete;

    /**
     * Runs `work` once the device has done all that comes before `after` (e.g. an event recorded
     * after the work that writes what `work` reads), and after the work posted before it. Without
     * an event, the work counts as ready right away.
     */
    void post(Event after, std::function<void()> work);

    /**
     * Room for the data of the next piece of work, which the device writes on the stream of the
     * runtime (e.g. with gathers): a buffer of pinned host memory of its own while it runs; while a
     * graph gets recorded, device memory of the recording.
     */
    Staging take(StreamRuntime& runtime);

    /**
     * Runs `work` with the data the device has written into `staging` (see take()) by now on the
     * stream of the runtime, once it is there; while a graph gets recorded, after each launch of
     * it.
     */
    void post(StreamRuntime& runtime,
              const Staging& staging,
              std::function<void(const void* data)> work);

    /**
     * Hands over the work of this queue that is ready; returns whether there was any.
     */
    bool progress();

    /**
     * Waits until all work posted to this queue so far is done, for the device as long as needed.
     */
    void drain();

    /// whether all work posted so far is done
    [[nodiscard]] bool empty() const { return pending() == 0; }

    /// the work posted so far that is not done yet
    [[nodiscard]] std::size_t pending() const {
      return posted_ - done_.load(std::memory_order_acquire);
    }

private:
    struct Item {
      Event after;
      std::function<void()> work;
    };

    // the work that waits for the device
    std::deque<Item> items_;
    std::size_t posted_{0};
    // the work done so far, on whatever thread it ran
    std::atomic<std::size_t> done_{0};
    // where the ready work becomes actions (tasks with a task graph), one after the other
    HostRuntime runtime_;
    HostActor& actor_;
    std::size_t bytes_;
    std::unique_ptr<HostBufferPool> buffers_;
  };

  /**
   * The work posted to the queues while a graph gets recorded (e.g. a super-timestep, see
   * beginRecording()). The device writes the data of the work into device memory of the recording,
   * at addresses that the graph keeps (see stage()). After each launch of the graph, launched()
   * copies all of it to a buffer of pinned host memory of its own, in the order of the stream, and
   * posts the work to its queues, in the order it was recorded, for once the copy is there. So the
   * launches may follow each other right away.
   */
  class Recording {
public:
    Recording() = default;
    ~Recording();

    Recording(const Recording&) = delete;
    Recording& operator=(const Recording&) = delete;
    Recording(Recording&&) = delete;
    Recording& operator=(Recording&&) = delete;

    /**
     * Device memory for `bytes` of data. It may close the node that the runtime captures from its
     * stream (see StreamRuntime::pauseCapture()), as no memory may be allocated while a stream gets
     * captured.
     */
    void* stage(std::size_t bytes, StreamRuntime& runtime);

    /// records `work` of `queue`, which reads the data at `staged` (from stage())
    void add(Queue& queue, const void* staged, std::function<void(const void* data)> work);

    /**
     * After a launch of the graph on the stream of the runtime: copies the data to the host, and
     * posts the work.
     */
    void launched(StreamRuntime& runtime);

    /// whether any work got recorded
    [[nodiscard]] bool empty() const { return items_.empty(); }

    /// the device memory of the recording
    [[nodiscard]] std::size_t bytes() const;

private:
    struct Chunk {
      void* device{nullptr};
      std::size_t capacity{0};
      std::size_t used{0};
    };
    struct Item {
      Queue* queue{nullptr};
      std::size_t chunk{0};
      std::size_t offset{0};
      std::function<void(const void*)> work;
    };

    std::vector<Chunk> chunks_;
    std::vector<Item> items_;
    // where the chunks are in a copy, and the copies (shared with the work that reads them)
    std::vector<std::size_t> copyOffsets_;
    std::shared_ptr<HostBufferPool> copies_;
  };

  /**
   * An actor that asks the device whether an event has completed; or, without a device, takes all
   * events as completed.
   */
  HostActor();

  /// an actor that asks `completed` instead (e.g. for tests)
  explicit HostActor(Completed completed);

  ~HostActor() = default;

  HostActor(const HostActor&) = delete;
  HostActor& operator=(const HostActor&) = delete;
  HostActor(HostActor&&) = delete;
  HostActor& operator=(HostActor&&) = delete;

  /// the actor of the process
  static HostActor& instance();

  /**
   * Hands the ready work to `executor` from now on, which runs what it gets in that order on a
   * thread of its own (e.g. the thread of the outputs); an empty one lets the queues run it again.
   * Only while no work is pending.
   */
  void setExecutor(Executor executor);

  /// whether an executor takes the ready work
  [[nodiscard]] bool hasExecutor() const { return static_cast<bool>(executor_); }

  /**
   * Hands over the ready work of all queues; returns whether there was any. Does nothing while a
   * graph gets recorded.
   */
  bool progress();

  /**
   * Waits until all work posted so far is done, for the device as long as needed.
   */
  void drain();

  /// the work posted to all queues that is not done yet
  [[nodiscard]] std::size_t pending() const;

  /**
   * The work that the queues get from now on goes into the recording (e.g. while a super-timestep
   * gets recorded), until endRecording().
   */
  void beginRecording(Recording& recording);
  void endRecording();

  private:
  [[nodiscard]] bool ready(const Event& event) const;

  Completed completed_;
  Executor executor_;
  std::vector<Queue*> queues_;
  Recording* recording_{nullptr};
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_HOSTACTOR_H_
