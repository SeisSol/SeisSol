#ifndef INITIALIZER_BOUNDARY_H_
#define INITIALIZER_BOUNDARY_H_

#include "Initializer/typedefs.hpp"
#include "Initializer/tree/LTSTree.hpp"
#include "Parallel/Helper.hpp"

#ifndef ACL_DEVICE
# define MEMKIND_BOUNDARY  initializer::AllocationMode::HostOnly
#else
# define MEMKIND_BOUNDARY  useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit
#endif // ACL_DEVICE

namespace seissol {
  namespace initializer {
    struct Boundary;
  }
}

struct seissol::initializer::Boundary {
  Variable<BoundaryFaceInformation> faceInformation;
  
  void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(faceInformation, mask, 1, MEMKIND_BOUNDARY);
  }
};
#endif
