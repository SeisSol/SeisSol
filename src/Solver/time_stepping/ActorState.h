#ifndef SEISSOL_ACTORSTATE_H
#define SEISSOL_ACTORSTATE_H
#include <queue>
#include <variant>
#include <omp.h>
#include <memory>
#include <algorithm>

namespace seissol {
namespace time_stepping {


struct AdvancedPredictionTimeMessage {
  double time;
};
struct AdvancedCorrectionTimeMessage {
  double time;
};

using Message = std::variant<AdvancedPredictionTimeMessage, AdvancedCorrectionTimeMessage>;

class MessageQueue {
 private:
  std::queue<Message> queue;
  omp_lock_t lock;

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
    Message message;
    omp_set_lock(&lock);
    message = queue.front();
    queue.front();
    queue.pop();
    omp_unset_lock(&lock);
    return message;
  }

  bool hasMessages() {
    return !queue.empty();
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

  double nextCorrectionTime(double syncTime) const {
    return std::min(syncTime, correctionTime + maxTimeStepSize);
  }

  //! Returns time step s.t. we won't miss the sync point
  double timeStepSize(double syncTime) const {
    assert(correctionTime < syncTime);
    return std::min(syncTime - correctionTime, maxTimeStepSize);
  }
};

struct NeighborCluster {
  ClusterTimes ct;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;

  NeighborCluster(double maxTimeStepSize) {
    ct.maxTimeStepSize = maxTimeStepSize;
  }


};

template<class T> struct always_false : std::false_type {};
}
}

#endif //SEISSOL_ACTORSTATE_H
