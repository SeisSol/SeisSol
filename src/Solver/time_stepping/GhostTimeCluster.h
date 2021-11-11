#ifndef SEISSOL_GHOSTTIMECLUSTER_H
#define SEISSOL_GHOSTTIMECLUSTER_H

#include <list>
#include "Initializer/typedefs.hpp"
#include "AbstractTimeCluster.h"

namespace seissol::time_stepping {

class GhostTimeCluster : public AbstractTimeCluster {
 private:
  const int globalClusterId;
  const int otherGlobalClusterId;
  const MeshStructure* meshStructure;
  std::list<MPI_Request*> sendQueue;
  std::list<MPI_Request*> receiveQueue;

  void sendCopyLayer();
  void receiveGhostLayer();

  bool testQueue(std::list<MPI_Request*>& queue);
  bool testForCopyLayerSends();
  bool testForGhostLayerReceives();

 public:
  GhostTimeCluster(double maxTimeStepSize,
                   int timeStepRate,
                   double timeTolerance,
                   int globalTimeClusterId,
                   int otherGlobalTimeClusterId,
                   const MeshStructure* meshStructure
  );
  ActResult act() override;

  void start() override;
  void predict() override;
  void correct() override;
  bool mayPredict() override;
  bool mayCorrect() override;
  bool maySync() override;
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
  void reset() override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  double lastSendTime = -1.0;
  double lastReceiveTime = -1.0;

};


}



#endif //SEISSOL_GHOSTTIMECLUSTER_H
