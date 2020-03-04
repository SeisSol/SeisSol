#ifndef SEISSOL_ACTORSTATE_H
#define SEISSOL_ACTORSTATE_H
#include <queue>
#include <variant>
#include <omp.h>
#include <memory>

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
    Synced,
    Finished
  };


struct ClusterTimes {
  double predictionTime = 0.0;
  double correctionTime = 0.0;
  double timeStepSize = std::numeric_limits<double>::infinity();

  // TODO(Lukas) Sync points
  double nextCorrectionTime() const {
    return correctionTime + timeStepSize;
  }
};

struct NeighborCluster {
  ClusterTimes ct;
  std::shared_ptr<MessageQueue> inbox = nullptr;
  std::shared_ptr<MessageQueue> outbox = nullptr;

  NeighborCluster(double timeStepSize) {
    ct.timeStepSize = timeStepSize;
  }


};


}
}

#endif //SEISSOL_ACTORSTATE_H
