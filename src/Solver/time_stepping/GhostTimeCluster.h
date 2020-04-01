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

 public:
  GhostTimeCluster(double maxTimeStepSize,
      double timeTolerance,
      int globalTimeClusterId,
      const MeshStructure* meshStructure)
      : AbstractTimeCluster(maxTimeStepSize, timeTolerance),
        globalClusterId(globalTimeClusterId),
        meshStructure(meshStructure) {}
  bool act() override;

  void predict() override;
  void correct() override;
  bool mayPredict() override;
  bool mayCorrect() override;
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
};


}



#endif //SEISSOL_GHOSTTIMECLUSTER_H
