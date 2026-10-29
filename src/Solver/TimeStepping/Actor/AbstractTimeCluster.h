// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_ABSTRACTTIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_ABSTRACTTIMECLUSTER_H_

#include "Memory/Tree/Layer.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Actor/StepParams.h"

#include <array>
#include <chrono>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace seissol::parallel::runtime {
class StreamRuntime;
} // namespace seissol::parallel::runtime

namespace seissol::solver {

class AbstractTimeCluster {
  private:
  ActorPriority priority_ = ActorPriority::Low;
  std::chrono::steady_clock::time_point timeOfLastStageChange_;
  bool alreadyPrintedTimeOut_ = false;

  //! the progress this cluster shows to the clusters waiting for it
  ActorProgress progress_;

  protected:
  //! how long the cluster may go without changing its state before it reports (see trackProgress)
  std::chrono::steady_clock::duration progressTimeout_{std::chrono::minutes(15)};

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
   * With concurrent clusters, makes the next action wait for the latest actions of the neighbors:
   * with access events, for their latest writes if it reads shared data, and for their latest reads
   * if it writes shared data (see DataAccess); otherwise, for their latest actions of any kind.
   */
  void waitForNeighbors(DataAccess access);

  /**
   * With concurrent clusters, publishes the event of the action just enqueued.
   */
  void publishEvent(DataAccess access);

  /**
   * What the device work of an action does with the shared data; the cell layers write it when
   * they predict, and read it when they correct.
   */
  [[nodiscard]] virtual DataAccess dataAccess(ActorAction action) const {
    return action == ActorAction::Predict ? DataAccess::Write : DataAccess::Read;
  }

  /// dataAccess(), or reads and writes of everything without access events
  [[nodiscard]] DataAccess actionAccess(ActorAction action) const;

  /// the event that work reading the data of `neighbor` waits for (see waitForNeighbors())
  [[nodiscard]] ActorEvent writesOf(const NeighborCluster& neighbor) const;
  /// the event that work overwriting data that `neighbor` reads waits for
  [[nodiscard]] ActorEvent readsOf(const NeighborCluster& neighbor) const;

  /**
   * Returns an event that completes with the device work enqueued so far; empty if there is none.
   */
  virtual ActorEvent recordActionEvent() { return {}; }

  /**
   * Makes the work enqueued from now on wait for the event.
   */
  virtual void waitForEvent(const ActorEvent& /*event*/) {}

  [[nodiscard]] bool concurrent() const { return concurrent_; }

  /**
   * Called after the time of this cluster has been set from outside.
   */
  virtual void timeSet(double /*time*/) {}

  /**
   * The host part of an action: decides and keeps the books on the host, before predict() or
   * correct() enqueue the device part of the action. It must not depend on the results of any
   * device work.
   */
  virtual StepWork prepare(ActorAction /*action*/) { return {}; }

  ActorState state_ = ActorState::Synced;
  ClusterTimes ct_;
  std::vector<NeighborCluster> neighbors_;
  //! runs the work of the actions on the host
  parallel::runtime::HostRuntime hostRuntime_;
  //! the names of the predictions and the corrections in a trace, once named
  std::array<std::uint32_t, 2> traceNames_{};
  bool traceNamed_{false};
  double syncTime_ = 0.0;
  //! the synchronization point before syncTime_, or the time the cluster started with
  double lastSyncTime_ = 0.0;
  StepContext stepContext_;

  [[nodiscard]] double timeStepSize() const;

  /**
   * A time step that equals the cluster's own step size recurs for the whole run, so a compute
   * graph recorded for it pays off. The truncated step taken right before a synchronization
   * point does not recur: caching it would add one graph per synchronization point.
   */
  [[nodiscard]] bool isRecurringTimestep(double timestep) const {
    return timestep == ct_.maxTimeStepSize;
  }

  /**
   * Parameters of the step that starts at the current correction time. Valid between reset() and
   * the next synchronization point.
   */
  [[nodiscard]] StepParams stepParams() const;

  void unsafePerformAction(ActorAction action);
  AbstractTimeCluster(double maxTimeStepSize, long timeStepRate, Executor executor);

  /**
   * Where an action of the host begins, in a trace (see monitoring::Trace); negative if the trace
   * does not record. The first one names the spans of the cluster in the trace.
   */
  double traceBegin();

  /// records the action of the host that began at `begin` in a trace
  void traceEnd(ActorAction action, double begin);

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

  //! whether a timeout aborts the simulation, for the clusters that can fail (see timeoutFail())
  bool abortOnTimeout_{true};

  //! whether someone wants a report (see requestReport())
  std::atomic<bool> reportWanted_{false};

  virtual void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate);

  public:
  /**
   * Reports the state of the cluster and of the neighbors it knows of, e.g. for a stall; on the
   * thread the cluster acts on.
   */
  virtual void reportState();

  /**
   * Makes the cluster report its state when it next looks at its progress, on the thread it acts
   * on; from any thread.
   */
  void requestReport() { reportWanted_ = true; }

  /**
   * A value that changes with every action of the cluster (the sum of the times it has published
   * to its neighbors), for watching it from another thread.
   */
  [[nodiscard]] double progressValue() const {
    return progress_.predictionTime + progress_.correctionTime;
  }

  /**
   * The time of the cluster as the device writes it back (see ClusterClock), for watching it from
   * another thread; null if the cluster has none.
   */
  [[nodiscard]] virtual const double* deviceTime() const { return nullptr; }

  /**
   * The actions whose work on the host is done (which, in a task graph, comes after the action),
   * for watching it from another thread.
   */
  [[nodiscard]] std::size_t hostWorkDone() const { return hostRuntime_.completedActions(); }

  /**
   * How long the cluster may go without changing its state before it reports, and whether a
   * timeout aborts the clusters that can fail (see the watchdog).
   */
  void setProgressTimeout(std::chrono::steady_clock::duration timeout, bool abort) {
    progressTimeout_ = timeout;
    abortOnTimeout_ = abort;
  }

  protected:
  [[nodiscard]] bool hasDifferentExecutorNeighbor() const;

  long timeStepRate_;

  //! only enqueue the device work of an action, instead of waiting for it
  bool concurrent_{false};
  //! with concurrent clusters, wait for the reads and writes of the neighbors apart
  bool accessEvents_{false};

  //! what the host part of the latest action has decided
  StepWork stepWork_{};

  //! enqueue the device work of the actions
  bool deviceWork_{true};
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

  /**
   * Reports when this cluster has not changed its state for too long, like act() does. For a
   * cluster that is not acted on since it has no action to take: in a deadlock, that holds for
   * all clusters.
   */
  void checkProgress() { trackProgress(false); }

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
   * With concurrent clusters, makes the actions wait only for the latest reads or writes of the
   * shared data by the neighbors that they depend on (see waitForNeighbors()).
   */
  void setAccessEvents(bool accessEvents) { accessEvents_ = accessEvents; }

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

  //! @brief The time steps the cluster has taken (corrected) so far.
  [[nodiscard]] long numberOfTimeSteps() const { return numberOfTimeSteps_; }

  /**
   * What the host part of the latest action has decided.
   */
  [[nodiscard]] const StepWork& lastStepWork() const { return stepWork_; }

  /**
   * Whether output samples fall into the next `steps` steps of this cluster; without any side
   * effect, and with the same decisions the host parts of these steps take.
   */
  [[nodiscard]] virtual bool outputsAhead(long /*steps*/) const { return false; }

  /**
   * Whether the next `steps` steps of this cluster all take its maximum time step size.
   */
  [[nodiscard]] bool regularStepsAhead(long steps) const;

  /**
   * Whether this cluster computes on the host, or exchanges data with clusters that do.
   */
  [[nodiscard]] virtual bool hostWork() const { return false; }

  /**
   * Lets the outputs decide about their samples when the work of a step runs, instead of in the
   * host part of the step; they then do not make a step irregular.
   */
  virtual void setRunTimeOutputs(bool /*runTimeOutputs*/) {}

  /**
   * The stream runtimes that the device work of this cluster goes through: its own, and those of
   * its outputs that take their samples on streams of their own. A recording of the work of the
   * cluster has to comprise them all.
   */
  [[nodiscard]] virtual std::vector<parallel::runtime::StreamRuntime*> runtimes() { return {}; }

  /**
   * Without device work, an action only runs its host part and keeps the books; its device work
   * then has to come from elsewhere, e.g. a recording.
   */
  void setDeviceWork(bool deviceWork) { deviceWork_ = deviceWork; }

  /**
   * The event that completes with the work of the latest action; empty if there is none.
   */
  [[nodiscard]] ActorEvent latestEvent() const { return progress_.event(); }

  /**
   * The size of the next step of this cluster.
   */
  [[nodiscard]] double nextTimeStepSize() const { return ct_.timeStepSize(syncTime_); }

  /**
   * Publishes the given event as the one that completes with the work of the latest action.
   */
  void publishEvent(ActorEvent event, DataAccess access = DataAccess::ReadWrite) {
    progress_.publishEvent(std::move(event), access);
  }

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

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_ABSTRACTTIMECLUSTER_H_
