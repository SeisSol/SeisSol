// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_COMMUNICATIONMANAGER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_COMMUNICATIONMANAGER_H_

#include "Parallel/HelperThread.h"
#include "Parallel/Pin.h"
#include "Solver/TimeStepping/Halo/GhostCluster.h"
#include "Solver/TimeStepping/Halo/TaskHaloTransport.h"

#include <atomic>
#include <memory>
#include <thread>
#include <vector>

namespace seissol::solver {
class AbstractCommunicationManager {
  public:
  using GhostClustersT = std::vector<std::unique_ptr<GhostCluster>>;
  virtual void progression() = 0;
  [[nodiscard]] virtual bool checkIfFinished() const = 0;
  virtual void reset(double newSyncTime);

  /**
   * With a task graph on the host: called by the main thread once it has created all tasks of the
   * epoch, before it waits for them. Keeps the exchanges going until all have completed, unless a
   * thread of its own does.
   */
  virtual void awaitExchanges() {}

  virtual ~AbstractCommunicationManager() = default;

  GhostClustersT* getGhostClusters();

  protected:
  explicit AbstractCommunicationManager(GhostClustersT ghostClusters);
  bool poll();
  GhostClustersT ghostClusters_;
};

class SerialCommunicationManager : public AbstractCommunicationManager {
  public:
  explicit SerialCommunicationManager(GhostClustersT ghostClusters);
  void progression() override;
  [[nodiscard]] bool checkIfFinished() const override;
};

class ThreadedCommunicationManager : public AbstractCommunicationManager {
  public:
  ThreadedCommunicationManager(GhostClustersT ghostClusters, const parallel::Pinning* pinning);
  void progression() override;
  [[nodiscard]] bool checkIfFinished() const override;
  void reset(double newSyncTime) override;

  private:
  seissol::parallel::HelperThread helper_;
};

/**
 * The halo exchange in a task graph on the host (see TaskHaloTransport). The main thread acts for
 * the ghost clusters, so that the tasks of the exchanges come in line with the tasks of the
 * clusters; their exchanges get started and tested by a thread of their own, or else by the main
 * thread, between the tasks it takes.
 */
class TaskCommunicationManager : public AbstractCommunicationManager {
  public:
  TaskCommunicationManager(GhostClustersT ghostClusters,
                           TaskHaloEngine& engine,
                           const parallel::Pinning* pinning,
                           bool thread);
  void progression() override;
  [[nodiscard]] bool checkIfFinished() const override;
  void reset(double newSyncTime) override;
  void awaitExchanges() override;

  private:
  TaskHaloEngine& engine_;
  // the main thread has created all tasks of the epoch
  std::atomic<bool> allCreated_{false};
  std::unique_ptr<seissol::parallel::HelperThread> helper_;
};

} // end namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_COMMUNICATIONMANAGER_H_
