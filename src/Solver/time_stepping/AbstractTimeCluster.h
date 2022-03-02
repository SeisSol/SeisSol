#ifndef SEISSOL_ACTOR_H
#define SEISSOL_ACTOR_H

#include <vector>
#include <memory>
#include <chrono>
#include "ActorState.h"

namespace seissol::time_stepping {

class AbstractTimeCluster {
private:
  ActorPriority priority = ActorPriority::Low;
protected:
  ActorState state = ActorState::Synced;
  ClusterTimes ct;
  std::vector<NeighborCluster> neighbors;
  double syncTime = 0.0;

  [[nodiscard]] double timeStepSize() const;

  void unsafePerformAction(ActorAction action);
public:
  AbstractTimeCluster(double maxTimeStepSize, long timeStepRate);

  virtual ~AbstractTimeCluster() = default;

  virtual ActorAction getNextLegalAction();
  virtual ActResult act();
  virtual bool mayPredict();
  virtual bool mayCorrect();
  virtual bool maySync();
  virtual void start() = 0;
  virtual void predict() = 0;
  [[nodiscard]] ActorState getState() const;
  virtual void correct() = 0;
  virtual bool processMessages();
  virtual void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) = 0;
  virtual void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) = 0;
  virtual void reset();
  virtual void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) = 0;

  ///* Returns the priority of the cluster. Larger numbers indicate a higher priority.
  ///* Can be used e.g. to always update copy clusters before interior ones.
  [[nodiscard]] virtual ActorPriority getPriority() const;
  virtual void setPriority(ActorPriority priority);

  void connect(AbstractTimeCluster& other);
  void updateSyncTime(double newSyncTime);
  [[nodiscard]] bool synced() const;

  void setPredictionTime(double time);
  void setCorrectionTime(double time);

  long timeStepRate;
  //! number of time steps
  long numberOfTimeSteps;
  std::chrono::steady_clock::time_point lastStateChange;
  const std::chrono::seconds timeout = std::chrono::minutes(15);
  bool alreadyPrintedTimeOut = false;

};

}



#endif //SEISSOL_ACTOR_H
