#pragma once

#include "Solver/time_stepping/GenericGhostTimeCluster.h"
#ifdef REQUIRED_COMM_LAYERS_PREFETCH
#include "Solver/time_stepping/GhostTimeClusterWithPrefetch.h"
#endif // REQUIRED_COMM_LAYERS_PREFETCH

namespace seissol::time_stepping {
struct GhostTimeClusterType {
  public:
#ifdef REQUIRED_COMM_LAYERS_PREFETCH
  using value = GhostTimeClusterWithPrefetch;
#else
  using value = GenericGhostTimeCluster;
#endif // REQUIRED_COMM_LAYERS_PREFETCH
};
} // namespace seissol::time_stepping