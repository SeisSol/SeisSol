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
#include <cstddef>
#include <memory>
#include <mutex>
#include <queue>
#include <variant>

namespace seissol::time_stepping {

struct AdvancedPredictionTimeMessage {
  double time;
  long stepsSinceSync;
};

struct AdvancedCorrectionTimeMessage {
  double time;
  long stepsSinceSync;
};

using Message = std::variant<AdvancedPredictionTimeMessage, AdvancedCorrectionTimeMessage>;

inline std::ostream& operator<<(std::ostream& stream, const Message& message);

class MessageQueue {
  private:
  std::queue<Message> queue_;
  std::mutex mutex_;
  // the number of messages in the queue, which the receiver polls without taking the lock, while
  // the sender may push to the queue
  std::atomic<std::size_t> size_{0};

  public:
  MessageQueue() = default;
  ~MessageQueue() = default;

  void push(const Message& message);

  Message pop();

  [[nodiscard]] bool hasMessages() const;

  [[nodiscard]] size_t size() const;
};

enum class ActorState { Corrected, Predicted, Synced };

enum class ActorAction { Nothing, Correct, Predict, Sync, RestartAfterSync };

std::string actorStateToString(ActorState state);

struct ClusterTimes {
  /**
   * The tolerance, in ticks (steps of the smallest time cluster), below which a remainder of time
   * before a synchronization point does not count as a tick of its own: the step before ends at the
   * synchronization point instead. It matches the tolerance with which the simulator compares its
   * times (TimeManager::getTimeTolerance()), and covers the rounding in the times of the clusters.
   */
  static constexpr double TickTolerance = 1e-5;

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

  /**
   * Whether the step that starts `steps` ticks after the last synchronization point is the last
   * one before the next.
   */
  [[nodiscard]] bool isLastStep(long steps) const { return steps + timeStepRate >= stepsUntilSync; }

  /**
   * The size of the step that starts at `time`, `steps` ticks after the last synchronization point.
   * The step counts decide which step is the last one before the synchronization point; that one
   * ends exactly there, all others are the maximum time step size.
   */
  [[nodiscard]] double timeStepSize(double time, long steps, double syncTime) const;

  //! The size of the step that starts at the correction time.
  [[nodiscard]] double timeStepSize(double syncTime) const;

  /**
   * The number of ticks from one synchronization point to the next, rounded up, but within the
   * tolerance (TickTolerance). All clusters count between the same synchronization times, not from
   * their own times, which may differ from them in the last bits: so they all agree on the number.
   */
  [[nodiscard]] long computeStepsUntilSyncTime(double oldSyncTime, double newSyncTime) const;

  //  [[nodiscard]] double& getTimeStepSize();

  [[nodiscard]] double getTimeStepSize() const { return maxTimeStepSize; }

  void setTimeStepSize(double newTimeStepSize) { maxTimeStepSize = newTimeStepSize; }
};

struct NeighborCluster {
  Executor executor;
  ClusterTimes ct;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;

  NeighborCluster(double maxTimeStepSize, int timeStepRate, Executor executor);
};

class DynamicRuptureScheduler {
  long lastCorrectionStepsInterior_ = -1;
  long lastCorrectionStepsCopy_ = -1;
  long numberOfDynamicRuptureFaces_;
  double outputTimestep_{};

  public:
  DynamicRuptureScheduler(long numberOfDynamicRuptureFaces, double outputTimestep);

  [[nodiscard]] bool mayComputeInterior(long curCorrectionSteps) const;

  void setLastCorrectionStepsInterior(long steps);

  void setLastCorrectionStepsCopy(long steps);

  [[nodiscard]] bool hasDynamicRuptureFaces() const;

  [[nodiscard]] double getOutputTimestep() const;
};

struct ActResult {
  bool isStateChanged = false;
};

enum class ActorPriority { Low, High };

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTORSTATE_H_
