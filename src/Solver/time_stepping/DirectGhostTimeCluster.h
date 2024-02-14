#pragma once

#include <list>
#include "Initializer/typedefs.hpp"
#include "Solver/time_stepping/AbstractGhostTimeCluster.h"


namespace seissol::time_stepping {
class DirectGhostTimeCluster : public AbstractGhostTimeCluster {
protected:
  virtual void sendCopyLayer();
  virtual void receiveGhostLayer();
  virtual bool testForGhostLayerReceives();

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

