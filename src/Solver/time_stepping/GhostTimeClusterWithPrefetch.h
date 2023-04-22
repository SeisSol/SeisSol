#pragma once

#include "Solver/time_stepping/GenericGhostTimeCluster.h"
#include <device.h>


namespace seissol::time_stepping {
class GhostTimeClusterWithPrefetch : public GenericGhostTimeCluster {
  public:
  GhostTimeClusterWithPrefetch(double maxTimeStepSize,
                               int timeStepRate,
                               int globalTimeClusterId,
                               int otherGlobalTimeClusterId,
                               const MeshStructure* meshStructure);
  ~GhostTimeClusterWithPrefetch();

  void sendCopyLayer() override;
  void receiveGhostLayer() override;

  bool testReceiveQueue();

  std::list<int> prefetchCopyLayer();
  void prefetchGhostRegion(int region);
  //void copyGhostRegionToDevice(int region);

  bool testForGhostLayerReceives() override;

  private:
  std::vector<void*> prefetchCopyRegionsStreams{};
  std::vector<void*> prefetchGhostRegionsStreams{};

  enum ReceiveState {
    RequiresMpiTesting,
    RequiresPrefetchTesting,
    Ready
  };
  std::vector<ReceiveState> receiveRegionsStates{};

  device::DeviceInstance& device = device::DeviceInstance::getInstance();
};
} // namespace seissol::time_stepping
