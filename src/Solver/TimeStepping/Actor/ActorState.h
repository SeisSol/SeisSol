// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_ACTORSTATE_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_ACTORSTATE_H_

#include "Common/Executor.h"
#include "Parallel/Runtime/EventPool.h"

#include <atomic>
#include <limits>
#include <mutex>
#include <string>
#include <utility>

namespace seissol::solver {

/**
 * An event that completes with the device work of an action, as it gets handed from a cluster to
 * the ones that wait for it.
 *
 * While it is held, an event from an event pool (e.g. the one of a stream runtime) stays reserved:
 * the pool does not hand it out to be recorded anew (see parallel::runtime::EventRef). So whoever
 * keeps an event to wait for it later -- a halo exchange that starts once the data is there, or a
 * replay that has to come after the work before it -- keeps the event itself. An event owned
 * elsewhere (e.g. by the halo exchange or by the recorder of super-timesteps) is not reserved; its
 * owner keeps it until nobody waits for it any more.
 */
class ActorEvent {
  public:
  ActorEvent() = default;

  /// an event owned elsewhere
  explicit ActorEvent(void* event) : event_(event) {}

  /// an event from an event pool; reserved while held
  explicit ActorEvent(parallel::runtime::EventRef event)
      : event_(event.get()), reference_(std::move(event)) {}

  [[nodiscard]] void* get() const { return event_; }

  explicit operator bool() const { return event_ != nullptr; }

  private:
  void* event_{nullptr};
  parallel::runtime::EventRef reference_;
};

/**
 * The progress of a cluster since the last synchronization point, as seen by the other clusters.
 *
 * Only the cluster itself writes its progress, after each of its actions; its neighbors read it,
 * possibly from other threads. The step counters are written last with release semantics, so that
 * a neighbor which has read them also sees the data of the steps they count.
 */
struct ActorProgress {
  std::atomic<long> predictionsSinceLastSync{0};
  std::atomic<long> stepsSinceLastSync{0};
  std::atomic<long> stepsUntilSync{0};
  std::atomic<double> predictionTime{0.0};
  std::atomic<double> correctionTime{0.0};

  /**
   * Makes `event` the one that completes with the device work of the latest action. Before the
   * step counters of the action, so that a neighbor which has seen the action also sees its event
   * (or a later one).
   */
  void publishEvent(ActorEvent event) {
    const std::scoped_lock lock(eventMutex_);
    event_ = std::move(event);
  }

  /**
   * The event that completes with the device work of the latest action; empty without concurrent
   * clusters. The copy keeps the event reserved while it is held.
   */
  [[nodiscard]] ActorEvent event() const {
    const std::scoped_lock lock(eventMutex_);
    return event_;
  }

  private:
  mutable std::mutex eventMutex_;
  ActorEvent event_;
};

enum class ActorState { Corrected, Predicted, Synced };

enum class ActorAction { Nothing, Correct, Predict, Sync, RestartAfterSync };

std::string actorStateToString(ActorState state);

struct ClusterTimes {
  double predictionTime = 0.0;
  double correctionTime = 0.0;
  double maxTimeStepSize = std::numeric_limits<double>::infinity();
  long stepsUntilSync = 0;
  long stepsSinceLastSync = 0;
  long predictionsSinceLastSync = 0;
  long predictionsSinceStart = 0;
  long stepsSinceStart = 0;
  long timeStepRate = -1;

  [[nodiscard]] double nextCorrectionTime(double syncTime) const;

  [[nodiscard]] long nextCorrectionSteps() const;

  //! Returns time step s.t. we won't miss the sync point
  [[nodiscard]] double timeStepSize(double syncTime) const;

  [[nodiscard]] long computeStepsUntilSyncTime(double oldSyncTime, double newSyncTime) const;

  //  [[nodiscard]] double& getTimeStepSize();

  [[nodiscard]] double getTimeStepSize() const { return maxTimeStepSize; }

  void setTimeStepSize(double newTimeStepSize) { maxTimeStepSize = newTimeStepSize; }
};

/**
 * When the data that a cluster provides for one of its steps becomes available to its neighbors.
 * Cell clusters (and the ghost clusters that stand in for remote ones) provide their data with
 * the prediction; face clusters provide theirs with their correction.
 */
enum class DataReadiness { AfterPrediction, AfterCorrection };

/**
 * What a cluster knows about one of the clusters it waits for.
 */
struct NeighborCluster {
  Executor executor;

  /// The times of the neighbor; the progress is taken over from `progress` before each decision.
  ClusterTimes ct;

  /// The progress the neighbor publishes.
  const ActorProgress* progress{nullptr};

  /// When the data of the neighbor for a step becomes available.
  DataReadiness dataReadiness{DataReadiness::AfterPrediction};

  NeighborCluster(double maxTimeStepSize, int timeStepRate, Executor executor);
};

struct ActResult {
  bool isStateChanged = false;
};

enum class ActorPriority { Low, High };

/**
 * What the host part of an action has decided about its device part.
 */
struct StepWork {
  /// the device part takes output samples in addition to the regular work of the step
  bool outputs{false};

  /// the cluster computes on the host, or exchanges data with clusters that do
  bool hostWork{false};

  /// such a step cannot be replayed from a recording of regular steps
  [[nodiscard]] bool irregular() const { return outputs || hostWork; }
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_ACTORSTATE_H_
