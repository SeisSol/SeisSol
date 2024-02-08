#pragma once

#include <list>
#include "Initializer/typedefs.hpp"
#include "AbstractTimeCluster.h"

namespace seissol::time_stepping {
class AbstractGhostTimeCluster : public AbstractTimeCluster {
  protected:
  const int globalClusterId;
  const int otherGlobalClusterId;
  const MeshStructure* meshStructure;
  std::list<unsigned int> sendQueue;
  std::list<unsigned int> receiveQueue;

  double lastSendTime = -1.0;

  virtual void sendCopyLayer() = 0;
  virtual void receiveGhostLayer() = 0;

  bool testQueue(MPI_Request* requests, std::list<unsigned int>& regions);
  bool testForCopyLayerSends();
  virtual bool testForGhostLayerReceives() = 0;

  void start() override;
  void predict() override;
  void correct() override;
  bool mayPredict() override;
  bool mayCorrect() override;
  bool maySync() override;
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  public:
  AbstractGhostTimeCluster(double maxTimeStepSize,
                          int timeStepRate,
                          int globalTimeClusterId,
                          int otherGlobalTimeClusterId,
                          const MeshStructure* meshStructure);

  void reset() override;
  ActResult act() override;
};
} // namespace seissol::time_stepping
