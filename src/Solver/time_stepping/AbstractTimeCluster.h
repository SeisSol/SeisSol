#ifndef SEISSOL_ACTOR_H
#define SEISSOL_ACTOR_H

#include <vector>
#include <memory>
#include <chrono>
#include "ActorState.h"

namespace seissol::time_stepping {

class AbstractTimeCluster {
private:
    int priority = 0;
protected:
  ActorState state = ActorState::Synced;
  ClusterTimes ct;
  std::vector<NeighborCluster> neighbors;
  double syncTime = 0.0;
  double timeTolerance;

  [[nodiscard]] double timeStepSize() const;

public:
  // TODO(Lukas) Move a lot of these to protected or private
  AbstractTimeCluster(double maxTimeStepSize, double timeTolerance, long timeStepRate);
  virtual ~AbstractTimeCluster() = default;

  virtual bool act();
  virtual bool mayPredict();
  virtual bool mayCorrect();
  virtual bool maySync();
  virtual void start() = 0;
  virtual void predict() = 0;
  ActorState getState() const;
  virtual void correct() = 0;
  virtual bool processMessages();
  virtual void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) = 0;
  virtual void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) = 0;
  virtual void reset();
  virtual void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) = 0;
  ///* Returns the priority of the cluster. Larger numbers indicate a higher priority.
  ///* Can be used e.g. to always update copy clusters before interior ones.
  virtual int getPriority() const;
  virtual void setPriority(int priority);

  void connect(AbstractTimeCluster& other);
  void updateSyncTime(double newSyncTime);
  [[nodiscard]] bool synced() const;
  long timeStepRate;
  //! number of time steps
  long numberOfTimeSteps;
  std::chrono::steady_clock::time_point lastStateChange;
  const std::chrono::seconds timeout = std::chrono::seconds(120);
  bool alreadyPrintedTimeOut = false;

};

}



#endif //SEISSOL_ACTOR_H
