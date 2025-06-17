// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_COMMUNICATIONMANAGER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_COMMUNICATIONMANAGER_H_

#include "Parallel/Pin.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include <atomic>
#include <memory>
#include <thread>
#include <vector>

namespace seissol::time_stepping {
class AbstractCommunicationManager {
  public:
  using GhostClustersT = std::vector<std::unique_ptr<AbstractGhostTimeCluster>>;
  virtual void progression() = 0;
  [[nodiscard]] virtual bool checkIfFinished() const = 0;
  virtual void reset(double newSyncTime);

  virtual ~AbstractCommunicationManager() = default;

  GhostClustersT* getGhostClusters();

  protected:
  explicit AbstractCommunicationManager(GhostClustersT ghostClusters);
  bool poll();
  GhostClustersT ghostClusters;
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

  ~ThreadedCommunicationManager() override;

  private:
  std::thread thread;
  std::atomic<bool> shouldReset;
  std::atomic<bool> isFinished;
  const parallel::Pinning* pinning;
};

} // end namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_COMMUNICATIONMANAGER_H_
