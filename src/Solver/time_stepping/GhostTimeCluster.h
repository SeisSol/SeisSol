#ifndef SEISSOL_GHOSTTIMECLUSTER_H
#define SEISSOL_GHOSTTIMECLUSTER_H

#include <list>
#include "Initializer/typedefs.hpp"
#include "AbstractTimeCluster.h"

namespace seissol::time_stepping {

class GhostTimeCluster : public AbstractTimeCluster {
 private:
  const int globalClusterId;
  const MeshStructure* meshStructure;
  std::list<MPI_Request*> sendQueue;
  std::list<MPI_Request*> receiveQueue;

  void sendCopyLayer();
  void receiveGhostLayer();
  bool testForCopyLayerSends();
  bool testForGhostLayerReceives();
  bool processMessages();

 public:
  GhostTimeCluster(double maxTimeStepSize,
      double timeTolerance,
      int globalTimeClusterId,
      const MeshStructure* meshStructure)
      : AbstractTimeCluster(maxTimeStepSize, timeTolerance),
        globalClusterId(globalTimeClusterId),
        meshStructure(meshStructure) {}
  bool act();

};

}



#endif //SEISSOL_GHOSTTIMECLUSTER_H
