// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_ABSTRACTTIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_ABSTRACTTIMECLUSTER_H_

#include "ActorState.h"
#include "Memory/Tree/Layer.h"
#include "Solver/TimeStepping/StepParams.h"

#include <chrono>
#include <memory>
#include <string>
#include <vector>

namespace seissol::time_stepping {

class AbstractTimeCluster {
  private:
  ActorPriority priority_ = ActorPriority::Low;
  std::chrono::steady_clock::time_point timeOfLastStageChange_;
  const std::chrono::seconds timeout = std::chrono::minutes(15);
  bool alreadyPrintedTimeOut_ = false;

  //! the progress this cluster shows to the clusters waiting for it
  ActorProgress progress_;

  protected:
  /**
   * Makes the current times of this cluster visible to the clusters waiting for it.
   */
  void publishProgress();

  /**
   * Takes over the progress the neighbors have published. For each neighbor that has predicted or
   * corrected since, the matching handler is called.
   */
  void refreshNeighbors();

  /**
   * Records whether an action has made progress; reports when there has been none for too long.
   */
  void trackProgress(bool progressed);

  /**
   * With concurrent clusters, makes the next action wait for the latest actions of the neighbors.
   */
  void waitForNeighbors();

  /**
   * With concurrent clusters, publishes the event of the action just enqueued.
   */
  void publishEvent();

  /**
   * Returns an event that completes with the device work enqueued so far; null if there is none.
   */
  virtual void* recordActionEvent() { return nullptr; }

  /**
   * Makes the work enqueued from now on wait for the event.
   */
  virtual void waitForEvent(void* /*event*/) {}

  [[nodiscard]] bool concurrent() const { return concurrent_; }

  /**
   * Called after the time of this cluster has been set from outside.
   */
  virtual void timeSet(double /*time*/) {}

  ActorState state_ = ActorState::Synced;
  ClusterTimes ct_;
  std::vector<NeighborCluster> neighbors_;
  double syncTime_ = 0.0;
  StepContext stepContext_;

  [[nodiscard]] double timeStepSize() const;

  /**
   * Parameters of the step that starts at the current correction time. Valid between reset() and
   * the next synchronization point.
   */
  [[nodiscard]] StepParams stepParams() const;

  void unsafePerformAction(ActorAction action);
  AbstractTimeCluster(double maxTimeStepSize, long timeStepRate, Executor executor);

  virtual bool mayPredict();
  virtual bool mayCorrect();

  virtual bool maySync();
  virtual void start() = 0;
  virtual void predict() = 0;
  virtual void correct() = 0;

  /**
   * Called when `neighbor` has predicted since the last refresh; its times are updated already.
   */
  virtual void handleNeighborPrediction(const NeighborCluster& neighbor) = 0;

  /**
   * Called when `neighbor` has corrected since the last refresh; its times are updated already.
   */
  virtual void handleNeighborCorrection(const NeighborCluster& neighbor) = 0;

  [[nodiscard]] virtual bool timeoutFail() const;

  virtual void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate);

  bool hasDifferentExecutorNeighbor();

  long timeStepRate_;

  //! only enqueue the device work of an action, instead of waiting for it
  bool concurrent_{false};
  //! number of time steps
  long numberOfTimeSteps_{0};
  Executor executor_;

  public:
  virtual ~AbstractTimeCluster() = default;

  virtual void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {}

  [[nodiscard]] virtual std::string description() const { return ""; }

  [[nodiscard]] Executor getExecutor() const;

  virtual ActorAction getNextLegalAction();
  virtual ActResult act();
  virtual void finalize();

  ///* Returns the priority of the cluster. Larger numbers indicate a higher priority.
  ///* Can be used e.g. to always update copy clusters before interior ones.
  [[nodiscard]] virtual ActorPriority getPriority() const;
  virtual void setPriority(ActorPriority priority);

  /**
   * With concurrent clusters, an action only enqueues its device work: it first makes it wait for
   * the device work of the latest action of each neighbor, and publishes an event for its own.
   */
  void setConcurrent(bool concurrent) { concurrent_ = concurrent; }

  /**
   * Connects two clusters that wait for each other.
   */
  void connect(AbstractTimeCluster& other);

  /**
   * Makes this cluster wait for `other`, without `other` waiting for this cluster in turn.
   *
   * All clusters need to be reset before any of them acts again after a synchronization point.
   */
  void observe(AbstractTimeCluster& other);
  void setSyncTime(double newSyncTime);

  [[nodiscard]] ActorState getState() const;
  [[nodiscard]] bool synced() const;
  virtual void reset();

  virtual void setTime(double time);

  virtual void finishPhase();

  [[nodiscard]] long getTimeStepRate() const;

  /**
   * The number of steps of the smallest time cluster until the synchronization point; set by
   * reset().
   */
  [[nodiscard]] long getStepsUntilSync() const;

  /**
   * When the data this cluster provides for a step becomes available to its neighbors.
   */
  [[nodiscard]] virtual DataReadiness dataReadiness() const;

  [[nodiscard]] std::string identifier() const {
    return description() + "-" + std::to_string(ct_.timeStepRate);
  }

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

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ABSTRACTTIMECLUSTER_H_
