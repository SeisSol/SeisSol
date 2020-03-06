#ifndef SEISSOL_ACTOR_H
#define SEISSOL_ACTOR_H

#include <vector>
#include <memory>
#include "ActorState.h"

namespace seissol {
namespace time_stepping {

class AbstractTimeCluster {
protected:
  ActorState state = ActorState::Synced;
  ClusterTimes ct;
  std::vector<NeighborCluster> neighbors;
  double syncTime = 0.0;
  double timeTolerance;

public:
  AbstractTimeCluster(double timeTolerance)
    : timeTolerance(timeTolerance) {}
  virtual ~AbstractTimeCluster() {}

  virtual bool act() = 0;

  void connect(AbstractTimeCluster& other) {
    neighbors.emplace_back(other.ct.maxTimeStepSize);
    other.neighbors.emplace_back(ct.maxTimeStepSize);
    neighbors.back().inbox = std::make_shared<MessageQueue>();
    other.neighbors.back().inbox = std::make_shared<MessageQueue>();
    neighbors.back().outbox = other.neighbors.back().inbox;
    other.neighbors.back().outbox = neighbors.back().inbox;
  }

  void updateSyncTime(double newSyncTime) {
    assert(newSyncTime > syncTime);
    syncTime = newSyncTime;
  }

  bool synced() const {
    return state == ActorState::Synced;
  }
};

}
}



#endif //SEISSOL_ACTOR_H
