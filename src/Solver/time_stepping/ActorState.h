#ifndef SEISSOL_ACTORSTATE_H
#define SEISSOL_ACTORSTATE_H

#include <omp.h>
#include <ostream>
#include <queue>
#include <limits>
#include <memory>
#include <string>
#include <variant>

namespace seissol::time_stepping {

class OmpLockGuard {
  private:
    omp_lock_t* lock;
  public:
  OmpLockGuard(omp_lock_t* lock) : lock(lock) {
   omp_set_lock(lock);
   };
  ~OmpLockGuard() {
    omp_unset_lock(lock);
  }
};

enum class ActorPriority {
  Low,
  High
};

struct AdvancedPredictionTimeMessage {
  double time;
  long stepsSinceSync;
};

struct AdvancedCorrectionTimeMessage {
  double time;
  long stepsSinceSync;
};

using Message = std::variant<AdvancedPredictionTimeMessage, AdvancedCorrectionTimeMessage>;
// Helper for std::visit variant pattern
template<class T> struct always_false : std::false_type {};

inline std::ostream& operator<<(std::ostream& stream, const Message& message);

class MessageQueue {
 private:
  std::queue<Message> queue;
  omp_lock_t lock;

 public:
  MessageQueue() {
    omp_init_lock(&lock);
   };
  ~MessageQueue() {
    omp_destroy_lock(&lock);
  }

  void push(Message const& message);

  Message pop();

  [[nodiscard]] bool hasMessages() const;

  [[nodiscard]] size_t size() const;
};

enum class ActorState {
  Corrected,
  Predicted,
  Synced
};

enum class ActorAction {
  Nothing,
  Correct,
  Predict,
  Sync,
  RestartAfterSync
};

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

  [[nodiscard]] long computeStepsUntilSyncTime(double oldSyncTime,
                                               double newSyncTime) const;

};

struct NeighborCluster {
  ClusterTimes ct;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;

  NeighborCluster(double maxTimeStepSize, int timeStepRate);

};

class DynamicRuptureScheduler {
  long lastCorrectionStepsInterior = -1;
  long lastCorrectionStepsCopy = -1;
  long lastFaultOutput = -1;
  long numberOfDynamicRuptureFaces;

  omp_lock_t drLock;

public:
  explicit DynamicRuptureScheduler(long numberOfDynamicRuptureFaces);

  [[nodiscard]] bool mayComputeInterior(long curCorrectionSteps) const;

  [[nodiscard]] bool mayComputeFaultOutput(long curCorrectionSteps) const;

  void setLastCorrectionStepsInterior(long steps);

  void setLastCorrectionStepsCopy(long steps);

  void setLastFaultOutput(long steps);

  [[nodiscard]] bool hasDynamicRuptureFaces() const;

  void lock();

  void unlock();
};

struct ActResult {
  bool isStateChanged = false;
};


}

#endif //SEISSOL_ACTORSTATE_H
