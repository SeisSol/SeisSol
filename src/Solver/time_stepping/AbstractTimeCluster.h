#ifndef SEISSOL_ACTOR_H
#define SEISSOL_ACTOR_H

#include "ActorState.h"
#include <Parallel/Runtime/Stream.hpp>
#include <chrono>
#include <memory>
#include <queue>
#include <vector>

namespace seissol::time_stepping {

class AbstractTimeCluster {
private:
  ActorPriority priority = ActorPriority::Low;
  std::chrono::steady_clock::time_point timeOfLastStageChange;
  const std::chrono::seconds timeout = std::chrono::minutes(15);
  bool alreadyPrintedTimeOut = false;

protected:
  seissol::parallel::runtime::StreamRuntime streamRuntime;

  ActorState state = ActorState{StateType::Synchronized,ComputeStep::Correct};
  ClusterTimes ct;
  std::vector<NeighborCluster> neighbors;
  double syncTime = 0.0;
  Executor executor;

  [[nodiscard]] double timeStepSize() const;

  AbstractTimeCluster(double maxTimeStepSize, long timeStepRate, Executor executor);

  bool maySynchronize();
  virtual void start() = 0;
  void synchronize();
  void restart();
  bool allComputed(ComputeStep step);
  void preCompute(ComputeStep step);
  virtual void runCompute(ComputeStep step) = 0;
  virtual bool pollCompute(ComputeStep step) { return true; }
  bool advanceState();
  void postCompute(ComputeStep step);
  virtual ComputeStep lastStep() const { return ComputeStep::Correct; }
  ComputeStep nextStep(ComputeStep step) const {
    switch (step) {
    case ComputeStep::Predict:
      return ComputeStep::Interact;
    case ComputeStep::Interact:
      return ComputeStep::Correct;
    case ComputeStep::Correct:
      return ComputeStep::Predict;
    }
    throw;
  }
  // needed, so that we can re-use events from empty steps
  virtual bool emptyStep(ComputeStep step) const {return true;}
  virtual bool processMessages();
  virtual void handleAdvancedComputeTimeMessage(ComputeStep step, const NeighborCluster& neighborCluster) = 0;
  virtual void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) = 0;

  bool hasDifferentExecutorNeighbor();


  long timeStepRate;
  //! number of time steps
  long numberOfTimeSteps;

  std::queue<void*> events;
public:
  virtual ~AbstractTimeCluster();

  virtual std::string description() const { return ""; }

  Executor getExecutor() const;

  virtual ActResult act();
  virtual void finalize();

  ///* Returns the priority of the cluster. Larger numbers indicate a higher priority.
  ///* Can be used e.g. to always update copy clusters before interior ones.
  [[nodiscard]] virtual ActorPriority getPriority() const;
  virtual void setPriority(ActorPriority priority);

  void connect(AbstractTimeCluster& other);
  void setSyncTime(double newSyncTime);

  [[nodiscard]] ActorState getState() const;
  [[nodiscard]] bool synchronized() const;
  virtual void reset();

  void setTime(double time);

  long getTimeStepRate();

  /**
   * @brief Returns the time step size of the cluster.
   * @return the time step size of the cluster.
   */
  double getClusterTimes();
  /**
   * @brief Sets the time step size of the cluster.
   * @param newTimeStepSize
   */
  void setClusterTimes(double newTimeStepSize);

  /**
   * @brief Returns the neighbor clusters of the cluster.
   * @return the pointer to the vector of neighbor clusters.
   */
  std::vector<NeighborCluster>* getNeighborClusters();

};

class CellCluster : public AbstractTimeCluster {
protected:
 bool emptyStep(ComputeStep step) const override {return step == ComputeStep::Interact;}
 ~CellCluster() override;
  CellCluster(double maxTimeStepSize, long timeStepRate, Executor executor);
public:
};

class FaceCluster : public AbstractTimeCluster {
protected:
 bool emptyStep(ComputeStep step) const override {return step != ComputeStep::Interact;}
 ~FaceCluster() override;
  FaceCluster(double maxTimeStepSize, long timeStepRate, Executor executor);
public:
};

} // namespace seissol::time_stepping



#endif //SEISSOL_ACTOR_H
