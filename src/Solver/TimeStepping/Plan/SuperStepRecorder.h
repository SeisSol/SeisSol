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
 * The clusters enqueue their work on streams of their own; further streams (e.g. of the halo
 * exchange) can take part as well. A recording starts on a stream of the recorder, forks all
 * streams from it, and joins them back at its end. A replay waits for the latest work of all
 * clusters and streams before the super-timestep, launches the graph on the stream of the
 * recorder, and makes all clusters and streams wait for it. Meanwhile, the clusters only run the
 * host parts of their actions.
 *
 * Only available on devices that can record graphs.
 */
class SuperStepRecorder {
  public:
  /// what makes two super-timesteps do the same work: the step sizes of all clusters
  using Key = std::vector<double>;

  SuperStepRecorder();
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
   * Starts recording the work the clusters enqueue from now on. Returns the event that starts the
   * recording: all clusters and streams wait for it, and, once replayed, it comes after all work
   * before the recording. It stays reserved while the copy is held.
   */
  ActorEvent beginRecording(const std::vector<AbstractTimeCluster*>& clusters,
                            const std::vector<void*>& streams);

  /**
   * Ends the recording, and keeps it for the key. It still needs to be replayed once to run.
   */
  void endRecording(const Key& key,
                    const std::vector<AbstractTimeCluster*>& clusters,
                    const std::vector<void*>& streams);

  /**
   * Ends the recording and launches it once, as replay() would; it is not kept for the key of the
   * super-timestep, but released once it has run. Its device time gets added to `deviceTime` (see
   * setTimer()).
   */
  void endRecordingAndLaunch(const std::vector<AbstractTimeCluster*>& clusters,
                             const std::vector<void*>& streams,
                             double* deviceTime = nullptr);

  /**
   * Notes the latest work of the clusters, before they run the host parts of a super-timestep.
   */
  void beginReplay(const std::vector<AbstractTimeCluster*>& clusters,
                   const std::vector<void*>& streams);

  /**
   * Replays the recording for the key after the work noted last. Its device time gets added to
   * `deviceTime` (see setTimer()).
   */
  void replay(const Key& key,
              const std::vector<AbstractTimeCluster*>& clusters,
              const std::vector<void*>& streams,
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
   * Releases the device resources; needs to happen before the device is finalized, and after
   * everybody has let go of the events of the recorder.
   */
  void dispose();

  private:
  DeviceTimer* timer_{nullptr};

#ifdef ACL_DEVICE
  /// records an event on the stream, from a pool that does not hand it out again while it is held
  ActorEvent recordEvent(void* stream);

  /// joins all streams back into the stream of the recorder, and ends the recording
  void endCapture(const std::vector<AbstractTimeCluster*>& clusters,
                  const std::vector<void*>& streams);

  /// launches the graph after the work noted last, and makes all clusters and streams wait for it
  void launch(const device::DeviceGraphHandle& graph,
              const std::vector<AbstractTimeCluster*>& clusters,
              const std::vector<void*>& streams,
              double* deviceTime);

  // the stream the recordings get replayed on, and the pool of the events around them
  parallel::runtime::StreamRuntime runtime_;
  void* stream_{nullptr};
  // completes with the latest replay
  ActorEvent lastEvent_;
  // starts the latest recording
  ActorEvent fork_;
  // the latest work before the super-timestep, until the replay has been enqueued after it
  std::vector<ActorEvent> waitFor_;
  device::DeviceGraphHandle recording_;
  std::map<Key, device::DeviceGraphHandle> graphs_;
  // the recordings launched once, until the event after each has completed
  std::vector<std::pair<device::DeviceGraphHandle, ActorEvent>> launchedOnce_;
#endif
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_PLAN_SUPERSTEPRECORDER_H_
