#ifndef SEISSOL_ACTORSTATE_H
#define SEISSOL_ACTORSTATE_H

#include <cassert>
#include <cmath>
#include <queue>
#include <variant>
#include <omp.h>
#include <memory>
#include <algorithm>
#include <type_traits>


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
// Helper for std::visit variant pattern
template<class T> struct always_false : std::false_type {};

inline std::ostream& operator<<(std::ostream& stream, const Message& message) {
  std::visit([&stream](auto&& msg) {
    using T = std::decay_t<decltype(msg)>;
    if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
      stream << "AdvancedPredictionTimeMessage, t = " << msg.time;
    } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
      stream << "AdvancedCorrectionTimeMessage, t = " << msg.time;
    } else {
      static_assert(always_false<T>::value, "non-exhaustive visitor");
    }
  }, message);
  return stream;
}


class MessageQueue {
 private:
  std::queue<Message> queue;
  omp_lock_t lock{};

 public:
  MessageQueue() {
    omp_init_lock(&lock);
  }
  ~MessageQueue() {
    omp_destroy_lock(&lock);
  }

  void push(Message const& message) {
    omp_set_lock(&lock);
    queue.push(message);
    omp_unset_lock(&lock);
  }

  Message pop() {
    omp_set_lock(&lock);
    const Message message = queue.front();
    queue.pop();
    omp_unset_lock(&lock);
    return message;
  }

  [[nodiscard]] bool hasMessages() const {
    return !queue.empty();
  }

  [[nodiscard]] size_t size() const {
      return queue.size();
  }
};

  enum class ActorState {
      Corrected,
      Predicted,
      Synced
};

  inline std::string actorStateToString(ActorState state) {
    switch (state) {
      case ActorState::Corrected:
        return "Corrected";
      case ActorState::Predicted:
        return "Predicted";
      case ActorState::Synced:
        return "Synced";
      }
    throw;
  }

struct ClusterTimes {
  double predictionTime = 0.0;
  double correctionTime = 0.0;
  double maxTimeStepSize = std::numeric_limits<double>::infinity();
  long stepsUntilSync = 0;
  long stepsSinceLastSync = 0;
  long predictionsSinceLastSync = 0;
  long timeStepRate = -1;

  [[nodiscard]] double nextCorrectionTime(double syncTime) const {
    return std::min(syncTime, correctionTime + maxTimeStepSize);
  }

  [[nodiscard]] long nextCorrectionSteps() const {
     return std::min(stepsSinceLastSync + timeStepRate, stepsUntilSync);
  }

  //! Returns time step s.t. we won't miss the sync point
  [[nodiscard]] double timeStepSize(double syncTime) const {
    assert(correctionTime < syncTime);
    return std::min(syncTime - correctionTime, maxTimeStepSize);
  }

  [[nodiscard]] long computeStepsUntilSyncTime(double oldSyncTime,
          double newSyncTime) const {
      const double timeDiff = newSyncTime-oldSyncTime;
      return static_cast<long>(std::ceil(timeStepRate*timeDiff/maxTimeStepSize));
  }

};

struct NeighborCluster {
  ClusterTimes ct;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;

  explicit NeighborCluster(double maxTimeStepSize, int timeStepRate) {
    ct.maxTimeStepSize = maxTimeStepSize;
    ct.timeStepRate = timeStepRate;
  }

};

}

#endif //SEISSOL_ACTORSTATE_H
