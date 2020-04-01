#ifndef SEISSOL_ACTOR_H
#define SEISSOL_ACTOR_H

#include <vector>
#include <memory>
#include "ActorState.h"

namespace seissol::time_stepping {

class AbstractTimeCluster {
protected:
  ActorState state = ActorState::Synced;
  ClusterTimes ct;
  std::vector<NeighborCluster> neighbors;
  double syncTime = 0.0;
  double timeTolerance;

  [[nodiscard]] double timeStepSize() const;

public:
  // TODO(Lukas) Move a lot of these to protected or private
  AbstractTimeCluster(double maxTimeStepSize, double timeTolerance);
  virtual ~AbstractTimeCluster() = default;

  virtual bool act();
  virtual bool mayPredict();
  virtual bool mayCorrect();
  virtual void predict() = 0;
  virtual void correct() = 0;
  virtual bool processMessages();
  virtual void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) = 0;
  virtual void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) = 0;

  void connect(AbstractTimeCluster& other);
  void updateSyncTime(double newSyncTime);
  [[nodiscard]] bool synced() const;
};

}



#endif //SEISSOL_ACTOR_H
