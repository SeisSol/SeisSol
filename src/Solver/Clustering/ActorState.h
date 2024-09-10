#ifndef SEISSOL_ACTORSTATE_H
#define SEISSOL_ACTORSTATE_H

#include <memory>
#include <mutex>
#include <queue>
#include <unordered_map>

#include "Common/Executor.h"

namespace seissol::solver::clustering {

enum class ComputeStep { Predict, Interact, Correct };

struct Message {
  ComputeStep step;
  double time;
  long stepsSinceSync;
  void* completionEvent;
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

enum class StateType { Synchronized, ComputeStart, ComputeDone };

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

  [[nodiscard]] long computeStepsUntilSyncTime(double oldSyncTime, double newSyncTime) const;

  double getTimeStepSize() const { return maxTimeStepSize; }

  void setTimeStepSize(double newTimeStepSize) { maxTimeStepSize = newTimeStepSize; }
};

struct NeighborCluster {
  Executor executor;
  ClusterTimes ct;
  ComputeStep lastStep;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;
  std::unordered_map<ComputeStep, void*> events;

  NeighborCluster(double maxTimeStepSize, int timeStepRate, Executor executor);
};

struct ActResult {
  bool isStateChanged = false;
};

enum class ActorPriority { Low, High };

} // namespace seissol::solver::clustering

#endif // SEISSOL_ACTORSTATE_H
