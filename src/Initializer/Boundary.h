#ifndef INITIALIZER_BOUNDARY_H_
#define INITIALIZER_BOUNDARY_H_

#include "IO/Instance/Checkpoint/CheckpointManager.hpp"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Initializer/typedefs.hpp"
#include "Parallel/Helper.hpp"

namespace seissol::initializer {

inline auto allocationModeBoundary() {
#ifndef ACL_DEVICE
  return AllocationMode::HostOnly;
#else
  return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
#endif
}

struct Boundary {
  Variable<BoundaryFaceInformation> faceInformation;

  void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(faceInformation, mask, 1, allocationModeBoundary());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) {}
};

} // namespace seissol::initializer

#endif
