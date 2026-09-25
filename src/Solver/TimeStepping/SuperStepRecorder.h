// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_SUPERSTEPRECORDER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_SUPERSTEPRECORDER_H_

#include "Solver/TimeStepping/AbstractTimeCluster.h"

#include <cstddef>
#include <map>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::time_stepping {

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
  ~SuperStepRecorder();

  SuperStepRecorder(const SuperStepRecorder&) = delete;
  SuperStepRecorder(SuperStepRecorder&&) = delete;
  SuperStepRecorder& operator=(const SuperStepRecorder&) = delete;
  SuperStepRecorder& operator=(SuperStepRecorder&&) = delete;

  [[nodiscard]] static bool available();

  [[nodiscard]] bool has(const Key& key) const;

  /**
   * Starts recording the work the clusters enqueue from now on.
   */
  void beginRecording(const std::vector<AbstractTimeCluster*>& clusters,
                      const std::vector<void*>& streams);

  /**
   * Ends the recording, and keeps it for the key. It still needs to be replayed once to run.
   */
  void endRecording(const Key& key,
                    const std::vector<AbstractTimeCluster*>& clusters,
                    const std::vector<void*>& streams);

  /**
   * Notes the latest work of the clusters, before they run the host parts of a super-timestep.
   */
  void beginReplay(const std::vector<AbstractTimeCluster*>& clusters,
                   const std::vector<void*>& streams);

  /**
   * Replays the recording for the key after the work noted last.
   */
  void replay(const Key& key,
              const std::vector<AbstractTimeCluster*>& clusters,
              const std::vector<void*>& streams);

  /**
   * The event that completes with the latest replay.
   */
  [[nodiscard]] void* lastEvent() const;

  /**
   * Releases the device resources; needs to happen before the device is finalized.
   */
  void dispose();

  private:
#ifdef ACL_DEVICE
  void* nextEvent();

  void* stream_{nullptr};
  std::vector<void*> events_;
  std::size_t eventIndex_{0};
  void* lastEvent_{nullptr};
  std::vector<void*> waitFor_;
  device::DeviceGraphHandle recording_;
  std::map<Key, device::DeviceGraphHandle> graphs_;
#endif
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_SUPERSTEPRECORDER_H_
