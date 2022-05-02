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
  std::chrono::steady_clock::time_point timeOfLastStageChange;
  const std::chrono::seconds timeout = std::chrono::minutes(15);
  bool alreadyPrintedTimeOut = false;

protected:
  ActorState state = ActorState::Synced;
  ClusterTimes ct;
  std::vector<NeighborCluster> neighbors;
  double syncTime = 0.0;

  [[nodiscard]] double timeStepSize() const;

  void unsafePerformAction(ActorAction action);
  AbstractTimeCluster(double maxTimeStepSize, long timeStepRate);

  virtual bool mayPredict();
  virtual bool mayCorrect();
  virtual bool maySync();
  virtual void start() = 0;
  virtual void predict() = 0;
  virtual void correct() = 0;
  virtual bool processMessages();
  virtual void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) = 0;
  virtual void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) = 0;
  virtual void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) = 0;


  long timeStepRate;
  //! number of time steps
  long numberOfTimeSteps;

public:
  virtual ~AbstractTimeCluster() = default;

  virtual ActorAction getNextLegalAction();
  virtual ActResult act();

  ///* Returns the priority of the cluster. Larger numbers indicate a higher priority.
  ///* Can be used e.g. to always update copy clusters before interior ones.
  [[nodiscard]] virtual ActorPriority getPriority() const;
  virtual void setPriority(ActorPriority priority);

  void connect(AbstractTimeCluster& other);
  void setSyncTime(double newSyncTime);

  [[nodiscard]] ActorState getState() const;
  [[nodiscard]] bool synced() const;
  virtual void reset();

  void setPredictionTime(double time);
  void setCorrectionTime(double time);

  long getTimeStepRate();

};

}



#endif //SEISSOL_ACTOR_H
