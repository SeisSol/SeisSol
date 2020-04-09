#ifndef SEISSOL_ACTORSTATE_H
#define SEISSOL_ACTORSTATE_H

#include <cassert>
#include <queue>
#include <variant>
#include <omp.h>
#include <memory>
#include <algorithm>
#include <type_traits>


namespace seissol::time_stepping {


struct AdvancedPredictionTimeMessage {
  double time;
};
struct AdvancedCorrectionTimeMessage {
  double time;
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
    //omp_destroy_lock(&lock);
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

  bool hasMessages() {
    return !queue.empty();
  }

  void clear() {
    omp_set_lock(&lock);
    while(!queue.empty()) {
      std::visit([](auto&& msg) {
        using T = std::decay_t<decltype(msg)>;
        if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
          assert(("Should never cancel AdvancedPredictionTime Message!", false)) ;
        } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
          // nop
        } else {
          static_assert(always_false<T>::value, "non-exhaustive visitor");
        }
      }, queue.front());
      std::cout << "Cancelled " << queue.front() << std::endl;
      queue.pop();
    }
    //queue = {};
    omp_unset_lock(&lock);
  }
};

  enum class ActorState {
    Corrected,
    Predicted,
    Synced
  };


struct ClusterTimes {
  double predictionTime = 0.0;
  double correctionTime = 0.0;
  double maxTimeStepSize = std::numeric_limits<double>::infinity();

  [[nodiscard]] double nextCorrectionTime(double syncTime) const {
    return std::min(syncTime, correctionTime + maxTimeStepSize);
  }

  //! Returns time step s.t. we won't miss the sync point
  [[nodiscard]] double timeStepSize(double syncTime) const {
    assert(correctionTime < syncTime);
    return std::min(syncTime - correctionTime, maxTimeStepSize);
  }
};

struct NeighborCluster {
  ClusterTimes ct;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;

  explicit NeighborCluster(double maxTimeStepSize) {
    ct.maxTimeStepSize = maxTimeStepSize;
  }


};

}

#endif //SEISSOL_ACTORSTATE_H
