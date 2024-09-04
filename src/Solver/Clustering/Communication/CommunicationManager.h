#ifndef SEISSOL_COMMUNICATIONMANAGER_H
#define SEISSOL_COMMUNICATIONMANAGER_H

#include "NeighborCluster.hpp"
#include "Parallel/Pin.h"
#include <Parallel/HelperThread.hpp>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <atomic>
#include <memory>
#include <thread>
#include <vector>

namespace seissol::solver::clustering::communication {
class AbstractCommunicationManager {
  public:
  using GhostClustersT = std::vector<std::unique_ptr<NeighborCluster>>;
  virtual void progression() = 0;
  [[nodiscard]] virtual bool checkIfFinished() const = 0;
  virtual void reset();

  virtual ~AbstractCommunicationManager() = default;

  GhostClustersT& getGhostClusters();

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
  void reset() override;

  private:
  parallel::HelperThread helper;
};

} // namespace seissol::solver::clustering::communication

#endif // SEISSOL_COMMUNICATIONMANAGER_H
