// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_ACTORSTATE_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_ACTORSTATE_H_

#include <Parallel/Runtime/Stream.h>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <string>
#include <unordered_map>

#include "Common/Executor.h"

namespace seissol::solver::clustering {

constexpr double PriorityHighest = 1.0;
constexpr double PriorityNormal = 0.5;
constexpr double PriorityLowest = 0.0;

enum class ComputeStep { Predict, Communicate, Interact, Correct };

struct Message {
  ComputeStep step;
  double time;
  long stepsSinceSync;
  std::optional<parallel::runtime::EventT> completionEvent;
};

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

enum class StateType { Waiting, Synchronized, ComputeStart, ComputeDone };

struct ActorState {
  StateType type;
  ComputeStep step;
};

std::string actorStateToString(ActorState state);

struct ClusterTimes {
  std::unordered_map<ComputeStep, double> time;
  std::unordered_map<ComputeStep, long> computeSinceStart;
  std::unordered_map<ComputeStep, long> computeSinceLastSync;
  double maxTimeStepSize = std::numeric_limits<double>::infinity();
  long stepsUntilSync = 0;
  long stepsSinceStart = 0;
  long timeStepRate = -1;

  [[nodiscard]] double nextComputeTime(ComputeStep step, double syncTime) const;

  [[nodiscard]] long nextSteps() const;

  //! Returns time step s.t. we won't miss the sync point
  [[nodiscard]] double timeStepSize(double syncTime) const;

  [[nodiscard]] std::optional<double> speculativeTimeStepSize(double syncTime, int lookahead) const;

  [[nodiscard]] long computeStepsUntilSyncTime(double oldSyncTime, double newSyncTime) const;

  double getTimeStepSize() const { return maxTimeStepSize; }

  void setTimeStepSize(double newTimeStepSize) { maxTimeStepSize = newTimeStepSize; }
};

struct NeighborCluster {
  Executor executor;
  ClusterTimes ct;
  ComputeStep lastStep;
  std::string identifier;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;
  std::unordered_map<ComputeStep, std::optional<parallel::runtime::EventT>> events;
  bool dependent = true;

  NeighborCluster(double maxTimeStepSize, int timeStepRate, Executor executor);
};

struct ActResult {
  bool isStateChanged = false;
};

} // namespace seissol::solver::clustering

#endif // SEISSOL_SRC_SOLVER_CLUSTERING_ACTORSTATE_H_
