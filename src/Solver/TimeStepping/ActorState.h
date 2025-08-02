// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTORSTATE_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTORSTATE_H_

#include <memory>
#include <mutex>
#include <queue>
#include <variant>

#include <Common/Executor.h>

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
  std::queue<Message> queue;
  std::mutex mutex;

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

struct NeighborCluster {
  Executor executor;
  ClusterTimes ct;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;

  NeighborCluster(double maxTimeStepSize, int timeStepRate, Executor executor);
};

class DynamicRuptureScheduler {
  long lastCorrectionStepsInterior = -1;
  long lastCorrectionStepsCopy = -1;
  long lastFaultOutput = -1;
  long numberOfDynamicRuptureFaces;
  double outputTimestep{};

  public:
  DynamicRuptureScheduler(long numberOfDynamicRuptureFaces, double outputTimestep);

  [[nodiscard]] bool mayComputeInterior(long curCorrectionSteps) const;

  [[nodiscard]] bool mayComputeFaultOutput(long curCorrectionSteps) const;

  void setLastCorrectionStepsInterior(long steps);

  void setLastCorrectionStepsCopy(long steps);

  void setLastFaultOutput(long steps);

  [[nodiscard]] bool hasDynamicRuptureFaces() const;

  [[nodiscard]] double getOutputTimestep() const;
};

struct ActResult {
  bool isStateChanged = false;
};

enum class ActorPriority { Low, High };

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTORSTATE_H_
