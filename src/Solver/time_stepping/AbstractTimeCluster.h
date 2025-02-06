// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIME_STEPPING_ABSTRACTTIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIME_STEPPING_ABSTRACTTIMECLUSTER_H_

#include "ActorState.h"
#include <chrono>
#include <memory>
#include <vector>

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
  AbstractTimeCluster(double maxTimeStepSize, long timeStepRate, Executor executor);

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

  bool hasDifferentExecutorNeighbor();


  long timeStepRate;
  //! number of time steps
  long numberOfTimeSteps;
  Executor executor;

public:
  virtual ~AbstractTimeCluster() = default;

  Executor getExecutor() const;

  virtual ActorAction getNextLegalAction();
  virtual ActResult act();
  virtual void finalize();

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

}




#endif // SEISSOL_SRC_SOLVER_TIME_STEPPING_ABSTRACTTIMECLUSTER_H_

