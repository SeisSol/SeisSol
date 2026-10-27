// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_PLAN_SUPERSTEPRECORDER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_PLAN_SUPERSTEPRECORDER_H_

#include "Monitoring/DeviceTimer.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"

#include <cstddef>
#include <map>
#include <set>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::solver {

/**
 * Records the device work of whole super-timesteps (steps of the largest cluster) into graphs, and
 * replays them.
 *
 * The work goes through the stream runtimes of the clusters, of their outputs and of the halo
 * exchange. While a super-timestep gets recorded, they all build one graph, node by node (see
 * parallel::runtime::StreamRuntime::attach()); where the work of one of them waits for that of
 * another one, the graph gets an edge. A replay launches the graph on the stream of the recorder,
 * after the latest work of the runtimes before; the runtimes then defer their wait for the replay
 * to the work they get next (see parallel::runtime::StreamRuntime::deferWait()). So one replay
 * after another costs one launch and one event, and nothing else. Meanwhile, the clusters only run
 * the host parts of their actions.
 *
 * Only available on devices that can build graphs node by node.
 */
class SuperStepRecorder {
  public:
  /// what makes two super-timesteps do the same work: the step sizes of all clusters
  using Key = std::vector<double>;

  using Runtimes = std::vector<parallel::runtime::StreamRuntime*>;

  /**
   * With `branches`, the independent branches of the work of the runtimes become nodes of their own
   * (see parallel::runtime::StreamRuntime::attach()).
   */
  explicit SuperStepRecorder(bool branches = true);
#ifdef ACL_DEVICE
  ~SuperStepRecorder();
#else
  ~SuperStepRecorder() = default;
#endif

  SuperStepRecorder(const SuperStepRecorder&) = delete;
  SuperStepRecorder(SuperStepRecorder&&) = delete;
  SuperStepRecorder& operator=(const SuperStepRecorder&) = delete;
  SuperStepRecorder& operator=(SuperStepRecorder&&) = delete;

  [[nodiscard]] static bool available();

  [[nodiscard]] bool has(const Key& key) const;

  /**
   * Starts recording the work the runtimes enqueue from now on. Inside the recording, the work
   * must not wait for anything from before it, which its launch waits for; so the clusters publish
   * no event until they have work in it.
   */
  void beginRecording(const std::vector<AbstractTimeCluster*>& clusters, const Runtimes& runtimes);

  /**
   * Ends the recording, and keeps it for the key. It still needs to be replayed once to run.
   */
  void endRecording(const Key& key);

  /**
   * Ends the recording and launches it once, as replay() would; it is not kept for the key of the
   * super-timestep, but released once it has run. Its device time gets added to `deviceTime` (see
   * setTimer()).
   */
  void endRecordingAndLaunch(const std::vector<AbstractTimeCluster*>& clusters,
                             const Runtimes& runtimes,
                             double* deviceTime = nullptr);

  /**
   * Notes the latest work of the runtimes, before the clusters run the host parts of a
   * super-timestep.
   */
  void beginReplay(const Runtimes& runtimes);

  /**
   * Replays the recording for the key after the work noted last. Its device time gets added to
   * `deviceTime` (see setTimer()).
   */
  void replay(const Key& key,
              const std::vector<AbstractTimeCluster*>& clusters,
              const Runtimes& runtimes,
              double* deviceTime = nullptr);

  /**
   * Times the replays with `timer`, from where they can start to where they have run.
   */
  void setTimer(DeviceTimer* timer) { timer_ = timer; }

  /**
   * The event that completes with the latest replay. It stays reserved while the copy is held.
   */
  [[nodiscard]] ActorEvent lastEvent() const;

  /**
   * Releases the recordings that have been neither recorded nor replayed since the last call (e.g.
   * at the end of each epoch), after the device is past them: super-timesteps whose step sizes
   * have changed would otherwise keep their recordings, and the memory of the device for them, up
   * to the end of the run.
   */
  void forgetUnused();

  /**
   * Releases the device resources; needs to happen before the device is finalized, and after
   * everybody has let go of the events of the recorder (the runtimes their deferred waits).
   */
  void dispose();

  private:
  DeviceTimer* timer_{nullptr};

#ifdef ACL_DEVICE
  /// launches the graph after the work noted last, and orders the work after it
  void launch(const device::DeviceGraphHandle& graph,
              const std::vector<AbstractTimeCluster*>& clusters,
              const Runtimes& runtimes,
              double* deviceTime);

  // builds the recordings, and launches them on its stream; the pool of the events around them
  parallel::runtime::StreamRuntime runtime_;
  void* stream_{nullptr};
  bool branches_{true};
  // completes with the latest replay
  ActorEvent lastEvent_;
  // the latest work before the super-timestep, until the replay has been enqueued after it
  std::vector<ActorEvent> waitFor_;
  std::map<Key, device::DeviceGraphHandle> graphs_;
  // the keys recorded or replayed since the last forgetUnused()
  std::set<Key> used_;
  // the recordings launched once, until the event after each has completed
  std::vector<std::pair<device::DeviceGraphHandle, ActorEvent>> launchedOnce_;
#endif
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_PLAN_SUPERSTEPRECORDER_H_
