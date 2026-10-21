// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_HOSTRUNTIME_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_HOSTRUNTIME_H_

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// tasks with dependencies on a list of runtime length (iterators), and detached tasks: OpenMP 5.0
#if defined(_OPENMP) && _OPENMP >= 201811
#define SEISSOL_HOST_TASK_GRAPH
#endif

namespace seissol::parallel::runtime {

/**
 * What the work of an action does with the data of its cluster, and with the data of the clusters
 * it observes (its neighbors). Besides, the work changes what only the cluster itself works on
 * (e.g. its DOFs).
 */
struct HostAccess {
  /// the work reads what the neighbors share with the cluster (e.g. their time data)
  bool readsNeighbors{false};

  /// the work writes what the cluster shares with its neighbors; it reads it in any case
  bool writesShared{false};
};

/**
 * How the clusters do their work on the host.
 */
enum class HostScheduling {
  /// The work of an action runs when the action is decided. A loop is a parallel region of its
  /// own, which splits the range among all threads evenly, as a static schedule does.
  Static,

  /// The main thread decides the actions of an epoch in one parallel region, and the work of each
  /// action becomes a task, which runs once the work it depends on is done; meanwhile, the main
  /// thread goes on deciding. A loop becomes tasks of a fixed number of elements each.
  TaskGraph
};

/**
 * The work of a cluster on the host: the counterpart of its StreamRuntime for the work of its
 * actions, and for the loops over its cells or faces therein.
 *
 * With a task graph, the work of an action may run long after the action has been decided, while
 * the main thread decides the next ones: so it has to take along by value what these change (e.g.
 * the planned output samples). The order of the work follows from what the work of each action
 * does with the data of its cluster and its neighbors (see HostAccess): the work that writes data
 * comes after all work that reads or writes it before, and the work that reads data after the work
 * that writes it before. Each runtime has two addresses that stand for its data in the
 * dependencies of the tasks: one for what only its cluster works on, one for what it shares.
 */
class HostRuntime {
  public:
  HostRuntime() = default;
  ~HostRuntime() = default;

  // the neighbors and the tasks refer to the runtime by the addresses of its members
  HostRuntime(const HostRuntime&) = delete;
  HostRuntime(HostRuntime&&) = delete;
  HostRuntime& operator=(const HostRuntime&) = delete;
  HostRuntime& operator=(HostRuntime&&) = delete;

  /// whether this build can run a task graph
  static constexpr bool taskGraphAvailable() {
#ifdef SEISSOL_HOST_TASK_GRAPH
    return true;
#else
    return false;
#endif
  }

  /**
   * Whether the OpenMP runtime goes on when it runs a detached task right away (e.g. with a single
   * thread, or with many tasks waiting), as LLVM's and Intel's runtimes do. Others wait for the
   * event of the task there; then, the thread that creates the tasks must not be the one that
   * fulfills their events.
   */
  static constexpr bool detachedTasksNeverBlock() {
#ifdef KMP_VERSION_MAJOR
    return true;
#else
    return false;
#endif
  }

  /**
   * Sets how all host runtimes schedule, and for a task graph the elements of a loop that each of
   * its tasks takes. Before any work is run.
   */
  static void setScheduling(HostScheduling scheduling, std::size_t taskSize) {
    settings().scheduling = taskGraphAvailable() ? scheduling : HostScheduling::Static;
    settings().taskSize = std::max(taskSize, std::size_t{1});
  }

  [[nodiscard]] static HostScheduling scheduling() { return settings().scheduling; }

  [[nodiscard]] static std::size_t taskSize() { return settings().taskSize; }

  /**
   * Sets how many actions of all runtimes may wait for their work in a task graph at most, before
   * throttle() holds back whoever decides the next ones; 0 for no limit.
   */
  static void setActionLimit(std::size_t limit) { settings().actionLimit = limit; }

  [[nodiscard]] static std::size_t actionLimit() { return settings().actionLimit; }

  /**
   * In a task graph, waits until fewer actions than the limit wait for their work; meanwhile runs
   * `progress()` (e.g. what keeps the halo exchange going) and takes tasks. So that the thread that
   * decides the actions stays not too far ahead of their work.
   */
  template <typename F>
  static void throttle(const F& progress) {
    const auto limit = actionLimit();
    while (limit > 0 && inTaskGraph() && outstandingActions() >= limit) {
      std::invoke(progress);
#ifdef SEISSOL_HOST_TASK_GRAPH
      if constexpr (detachedTasksNeverBlock()) {
        // the thread may need to keep the halo exchange going, and takes a task meanwhile
        yield();
      } else {
        // a task yield may do nothing at all (as with libgomp); waiting for the oldest work takes
        // tasks until that is done, and a thread of its own keeps the halo exchange going
        auto& oldest = settings().oldest;
        if (!oldest.empty()) {
          const char* token = oldest.front();
          oldest.pop_front();
#pragma omp taskwait depend(in : token[0])
        }
      }
#endif
    }
  }

  /**
   * The work of this cluster reads what `other` shares (if it says so in its access).
   */
  void observe(const HostRuntime& other) { neighbors_.push_back(&other.shared_); }

  /**
   * The address that stands for what the cluster shares, in the dependencies of the tasks; for
   * those that work on it besides the actions (e.g. the halo exchange).
   */
  [[nodiscard]] char* sharedToken() { return &shared_; }

  /**
   * The address that stands for what the neighbor `index` (in the order of observe()) shares.
   */
  [[nodiscard]] const char* neighborToken(std::size_t index) const { return neighbors_.at(index); }

  /**
   * Runs the work of an action: right away, or with a task graph as a task after the work before it
   * that it depends on by `access`. The work is taken by value.
   */
  template <typename F>
  void action(const HostAccess& access, F&& work) {
    created_.fetch_add(1, std::memory_order_relaxed);
#ifdef SEISSOL_HOST_TASK_GRAPH
    if (inTaskGraph()) {
      auto task = std::decay_t<F>(std::forward<F>(work));
      const char* const* neighbors = neighbors_.data();
      const std::size_t count = access.readsNeighbors ? neighbors_.size() : 0;
      outstanding().fetch_add(1, std::memory_order_relaxed);
      if constexpr (!detachedTasksNeverBlock()) {
        // for throttle(); the work of the cluster so far is done once the own token is free
        settings().oldest.push_back(&own_);
      }
      if (access.writesShared) {
#pragma omp task default(none) firstprivate(task)                                                  \
    depend(iterator(std::size_t i = 0 : count), in : neighbors[i][0])                              \
    depend(inout : this->own_, this->shared_)
        {
          std::invoke(task);
          done();
        }
      } else {
#pragma omp task default(none) firstprivate(task)                                                  \
    depend(iterator(std::size_t i = 0 : count), in : neighbors[i][0]) depend(in : this->shared_)   \
    depend(inout : this->own_)
        {
          std::invoke(task);
          done();
        }
      }
      return;
    }
#endif
    std::invoke(std::forward<F>(work));
    completed_.fetch_add(1, std::memory_order_release);
  }

  /// the actions whose work has been handed over so far; from any thread
  [[nodiscard]] std::size_t createdActions() const {
    return created_.load(std::memory_order_relaxed);
  }

  /// the actions whose work is done so far; from any thread
  [[nodiscard]] std::size_t completedActions() const {
    return completed_.load(std::memory_order_acquire);
  }

  /// the actions of all runtimes whose work runs as a task, and has not been done yet
  [[nodiscard]] static std::size_t outstandingActions() {
    return outstanding().load(std::memory_order_acquire);
  }

  /**
   * Lets the calling thread take a task meanwhile, if there is one to take (e.g. while the main
   * thread waits for the work to catch up with it).
   */
  static void yield() {
#ifdef SEISSOL_HOST_TASK_GRAPH
    if (inTaskGraph()) {
#pragma omp taskyield
    }
#endif
  }

  /**
   * Runs `body(begin, end)` on ranges that together cover [0, count) once each: one range per
   * thread, as a static schedule splits it; or, within the work of a task graph, one per task, of
   * taskSize() elements each. Whatever a range needs (e.g. scratch memory) can be set up in the
   * body, once per range. Returns once all ranges are done.
   */
  template <typename F>
  static void forRange(std::size_t count, const F& body) {
    if (count == 0) {
      return;
    }
#ifdef SEISSOL_HOST_TASK_GRAPH
    if (inTaskGraph()) {
      const std::size_t size = taskSize();
      const std::size_t tasks = (count + size - 1) / size;
#pragma omp taskloop default(none) shared(body) firstprivate(count, size, tasks) grainsize(1)
      for (std::size_t task = 0; task < tasks; ++task) {
        const auto begin = task * size;
        timed(body, begin, std::min(count, begin + size));
      }
      return;
    }
#endif
#ifdef _OPENMP
#pragma omp parallel default(none) shared(count, body)
    {
      const auto [begin, end] = staticRange(count,
                                            static_cast<std::size_t>(omp_get_thread_num()),
                                            static_cast<std::size_t>(omp_get_num_threads()));
      if (begin < end) {
        timed(body, begin, end);
      }
    }
#else
    timed(body, std::size_t{0}, count);
#endif
  }

  /**
   * The time all threads have spent in the ranges of the loops since resetLoopTime(), in seconds;
   * e.g. to tell how busy the threads were in an epoch.
   */
  [[nodiscard]] static double loopTime() {
    return static_cast<double>(loopNanoseconds().load(std::memory_order_relaxed)) * 1e-9;
  }

  static void resetLoopTime() { loopNanoseconds().store(0, std::memory_order_relaxed); }

  /**
   * Runs `body(i)` for every i in [0, count), split like forRange().
   */
  template <typename F>
  static void forEach(std::size_t count, const F& body) {
    forRange(count, [&body](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        std::invoke(body, i);
      }
    });
  }

  /**
   * Runs `body`, in which the main thread decides the actions of an epoch. With a task graph, it
   * runs in a parallel region whose other threads take the tasks meanwhile; it returns once all of
   * them are done.
   */
  template <typename F>
  static void epoch(F&& body) {
#ifdef SEISSOL_HOST_TASK_GRAPH
    if (scheduling() == HostScheduling::TaskGraph) {
#pragma omp parallel default(none) shared(body)
      {
#pragma omp masked
        {
          std::invoke(body);
#pragma omp taskwait
        }
      }
      settings().oldest.clear();
      return;
    }
#endif
    std::invoke(std::forward<F>(body));
  }

  /**
   * The range of thread `thread` of `threads` when [0, count) is split among them evenly: the
   * first `count % threads` threads take one element more than the others (as the static schedule
   * of LLVM's OpenMP runtime does).
   */
  static std::pair<std::size_t, std::size_t>
      staticRange(std::size_t count, std::size_t thread, std::size_t threads) {
    const auto base = count / threads;
    const auto extra = count % threads;
    const auto begin = thread * base + std::min(thread, extra);
    return {begin, begin + base + (thread < extra ? 1 : 0)};
  }

  private:
  /// whether the work runs as tasks: in the parallel region of an epoch with a task graph, also
  /// with a single thread (whose region does not count as active)
  [[nodiscard]] static bool inTaskGraph() {
#ifdef SEISSOL_HOST_TASK_GRAPH
    return scheduling() == HostScheduling::TaskGraph && omp_get_level() > 0;
#else
    return false;
#endif
  }

  /// runs a range of a loop, and counts its time
  template <typename F>
  static void timed(const F& body, std::size_t begin, std::size_t end) {
    const auto start = std::chrono::steady_clock::now();
    std::invoke(body, begin, end);
    const auto time = std::chrono::steady_clock::now() - start;
    loopNanoseconds().fetch_add(
        static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(time).count()),
        std::memory_order_relaxed);
  }

  static std::atomic<std::uint64_t>& loopNanoseconds() {
    static std::atomic<std::uint64_t> nanoseconds{0};
    return nanoseconds;
  }

  /// the work of an action that ran as a task is done
  void done() {
    completed_.fetch_add(1, std::memory_order_release);
    outstanding().fetch_sub(1, std::memory_order_release);
  }

  static std::atomic<std::size_t>& outstanding() {
    static std::atomic<std::size_t> outstanding{0};
    return outstanding;
  }

  struct Settings {
    HostScheduling scheduling{HostScheduling::Static};
    std::size_t taskSize{32};
    std::size_t actionLimit{0};
    // where a task yield may do nothing: the own tokens of the actions run as tasks, oldest first
    std::deque<const char*> oldest;
  };

  /// the settings of all host runtimes
  static Settings& settings() {
    static Settings settings;
    return settings;
  }

  // the addresses that stand for the data of the cluster in the dependencies of the tasks: what
  // only it works on, and what it shares; and those of what its neighbors share
  char own_{0};
  char shared_{0};
  std::vector<const char*> neighbors_;

  std::atomic<std::size_t> created_{0};
  std::atomic<std::size_t> completed_{0};
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_HOSTRUNTIME_H_
