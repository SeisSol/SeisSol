#ifndef VISCOELASTIC_INTEGRATIONDATA_H_
#define VISCOELASTIC_INTEGRATIONDATA_H_

#include "generated_code/tensor.h"

namespace seissol {
  namespace model {

    struct ViscoElasticLocalData {
      real sourceMatrix[seissol::tensor::ET::size()];
    };
    struct ViscoElasticNeighborData {
    };
  }
}

#endif
