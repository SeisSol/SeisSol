#pragma once

#include "Initializer/typedefs.hpp"
#include "Solver/Clustering/Communication/AbstractGhostTimeCluster.h"
#include <list>

namespace seissol::time_stepping {
class DirectGhostTimeCluster : public AbstractGhostTimeCluster {
  protected:
  void sendCopyLayer() override;
  void receiveGhostLayer() override;
  bool testForGhostLayerReceives() override;

  public:
  DirectGhostTimeCluster(double maxTimeStepSize,
                         int timeStepRate,
                         int globalTimeClusterId,
                         int otherGlobalTimeClusterId,
                         const MeshStructure* meshStructure,
                         bool persistent);
  void finalize() override;

  private:
  bool persistent;
};
} // namespace seissol::time_stepping
