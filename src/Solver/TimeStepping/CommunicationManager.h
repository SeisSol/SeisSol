// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_COMMUNICATIONMANAGER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_COMMUNICATIONMANAGER_H_

#include "Parallel/Pin.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include <Parallel/HelperThread.h>
#include <atomic>
#include <memory>
#include <thread>
#include <vector>

namespace seissol::time_stepping {
class AbstractCommunicationManager {
  public:
  using GhostClustersT = std::vector<std::unique_ptr<AbstractTimeCluster>>;
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

  private:
  seissol::parallel::HelperThread helper;
};

} // end namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_COMMUNICATIONMANAGER_H_
