#ifndef INITIALIZER_BOUNDARY_H_
#define INITIALIZER_BOUNDARY_H_

#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
// #include <generated_code/tensor.h>

namespace seissol {
  namespace initializers {
    struct Boundary;
  }
}

struct seissol::initializers::Boundary {
  Variable<BoundaryFaceInformation> faceInformation;
  
  void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(faceInformation, mask, 1, seissol::memory::Standard);
  }
};
#endif
