// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTORSTATE_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTORSTATE_H_

#include "Common/Executor.h"

#include <atomic>
#include <limits>
#include <string>

namespace seissol::time_stepping {

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

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTORSTATE_H_
