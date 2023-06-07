#ifndef INITIALIZER_BOUNDARY_H_
#define INITIALIZER_BOUNDARY_H_

#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>

#include "Initializer/tree/VariableContainer.hpp"


#ifndef ACL_DEVICE
# define MEMKIND_BOUNDARY  seissol::memory::Standard
#else
# define MEMKIND_BOUNDARY  seissol::memory::DeviceUnifiedMemory
#endif // ACL_DEVICE

namespace seissol {
  namespace initializers {
    struct Boundary;
  }
}

struct seissol::initializers::Boundary : seissol::initializers::LTSVariableContainer {
  Variable<BoundaryFaceInformation> faceInformation;
  
  virtual void addTo(LTSTree& tree) override {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(faceInformation, mask, 1, MEMKIND_BOUNDARY);
  }
};
#endif
