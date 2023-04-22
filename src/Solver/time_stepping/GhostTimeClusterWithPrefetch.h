#pragma once

#include "Solver/time_stepping/GenericGhostTimeCluster.h"

namespace seissol::time_stepping {
class GhostTimeClusterWithPrefetch : public GenericGhostTimeCluster {
  enum Direction {Send, Recv};
  public:
  GhostTimeClusterWithPrefetch(double maxTimeStepSize,
                               int timeStepRate,
                               int globalTimeClusterId,
                               int otherGlobalTimeClusterId,
                               const MeshStructure* meshStructure)
      : GenericGhostTimeCluster(maxTimeStepSize,
                                timeStepRate,
                                globalTimeClusterId,
                                otherGlobalTimeClusterId,
                                meshStructure) {}

  void sendCopyLayer() override;
  void receiveGhostLayer() override;
  bool testReceiveQueue();

  void syncCommRegion(int region, Direction type);

  bool testForGhostLayerReceives() override;
};
} // namespace seissol::time_stepping