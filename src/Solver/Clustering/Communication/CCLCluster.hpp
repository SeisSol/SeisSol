#pragma once

#include "Initializer/typedefs.hpp"
#include "Solver/Clustering/AbstractTimeCluster.h"
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/AbstractGhostTimeCluster.h>
#include <list>

namespace seissol::time_stepping {
class CCLCluster : public CellCluster {
  protected:
  const int globalClusterId;
  const int otherGlobalClusterId;

  double lastSendTime = -1.0;

  std::vector<RemoteCluster> copyClusters;
  std::vector<RemoteCluster> ghostClusters;

  virtual void sendCopyLayer();
  virtual void receiveGhostLayer();

  void start() override;
  void runCompute(ComputeStep step) override;
  void handleAdvancedComputeTimeMessage(ComputeStep step,
                                        const NeighborCluster& neighborCluster) override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  void* cclComm;

  public:
  CCLCluster(double maxTimeStepSize,
             int timeStepRate,
             int globalTimeClusterId,
             int otherGlobalTimeClusterId,
             const MeshStructure* meshStructure);

  void reset() override;
  ActResult act() override;
};
} // namespace seissol::time_stepping
