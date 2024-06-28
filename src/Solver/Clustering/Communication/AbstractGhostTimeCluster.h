#pragma once

#include "Initializer/typedefs.hpp"
#include "Solver/Clustering/AbstractTimeCluster.h"
#include <Solver/Clustering/ActorState.h>
#include <list>

namespace seissol::time_stepping {
struct RemoteCluster {
  void* data;
  std::size_t size;
  MPI_Datatype datatype;
  int rank;
  int tag;
};

class AbstractGhostTimeCluster : public CellCluster {
  protected:
  const int globalClusterId;
  const int otherGlobalClusterId;
  const MeshStructure* meshStructure;
  std::list<std::size_t> sendQueue;
  std::list<std::size_t> recvQueue;

  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> recvRequests;

  std::vector<RemoteCluster> copyClusters;
  std::vector<RemoteCluster> ghostClusters;

  double lastSendTime = -1.0;

  virtual void sendCopyLayer() = 0;
  virtual void receiveGhostLayer() = 0;

  bool testQueue(std::vector<MPI_Request>& requests, std::list<std::size_t>& remaining);
  bool testForCopyLayerSends();
  virtual bool testForGhostLayerReceives() = 0;

  void start() override;
  void runCompute(ComputeStep step) override;
  bool pollCompute(ComputeStep step) override;
  void handleAdvancedComputeTimeMessage(ComputeStep step,
                                        const NeighborCluster& neighborCluster) override;
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
