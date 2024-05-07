#ifndef MODEL_POROELASTIC_INTEGRATIONDATA_H_
#define MODEL_POROELASTIC_INTEGRATIONDATA_H_

#include "generated_code/tensor.h"

namespace seissol {
  namespace model {

    struct PoroelasticLocalData {
      real sourceMatrix[seissol::tensor::ET::size()];
      real G[NUMBER_OF_QUANTITIES];
      real typicalTimeStepWidth;
      real Zinv[NUMBER_OF_QUANTITIES][CONVERGENCE_ORDER*CONVERGENCE_ORDER];
    };
    struct PoroelasticNeighborData {};
  }
}

#endif
