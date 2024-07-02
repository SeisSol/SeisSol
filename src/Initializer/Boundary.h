#ifndef INITIALIZER_BOUNDARY_H_
#define INITIALIZER_BOUNDARY_H_

#include "IO/Instance/Checkpoint/CheckpointManager.hpp"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/typedefs.hpp"
#include "Parallel/Helper.hpp"

#ifndef ACL_DEVICE
# define MEMKIND_BOUNDARY  initializer::AllocationMode::HostOnly
#else
# define MEMKIND_BOUNDARY  useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit
#endif // ACL_DEVICE

namespace seissol::initializer {
  struct Boundary {
    Variable<BoundaryFaceInformation> faceInformation;
    
    void addTo(LTSTree& tree) {
      LayerMask mask = LayerMask(Ghost);
      tree.addVar(faceInformation, mask, 1, MEMKIND_BOUNDARY);
    }

    void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager, LTSTree* tree) {
    }
  };
} // namespace seissol::initializer
#endif
