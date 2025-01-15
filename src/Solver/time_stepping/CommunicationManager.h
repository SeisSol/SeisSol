// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIME_STEPPING_COMMUNICATIONMANAGER_H_
#define SEISSOL_SRC_SOLVER_TIME_STEPPING_COMMUNICATIONMANAGER_H_

#include <atomic>
#include <memory>
#include <thread>
#include <vector>
#include "Parallel/Pin.h"
#include "Solver/time_stepping/AbstractGhostTimeCluster.h"


namespace seissol::time_stepping {
class AbstractCommunicationManager {
public:
  using ghostClusters_t = std::vector<std::unique_ptr<AbstractGhostTimeCluster>>;
  virtual void progression() = 0;
  [[nodiscard]] virtual bool checkIfFinished() const = 0;
  virtual void reset(double newSyncTime);

  virtual ~AbstractCommunicationManager() = default;

  ghostClusters_t* getGhostClusters();

protected:
  explicit AbstractCommunicationManager(ghostClusters_t ghostClusters);
  bool poll();
  ghostClusters_t ghostClusters;

};

class SerialCommunicationManager : public AbstractCommunicationManager {
public:
  explicit SerialCommunicationManager(ghostClusters_t ghostClusters);
  void progression() override;
  [[nodiscard]] bool checkIfFinished() const override;
};

class ThreadedCommunicationManager : public AbstractCommunicationManager {
public:
  ThreadedCommunicationManager(ghostClusters_t ghostClusters,
                               const parallel::Pinning* pinning);
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


#endif // SEISSOL_SRC_SOLVER_TIME_STEPPING_COMMUNICATIONMANAGER_H_

