#ifndef MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_
#define MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_

#include <generated_code/tensor.h>

namespace seissol {
  namespace model {

    struct ViscoElasticLocalData {
      real E[tensor::E::size()];
      real w[tensor::w::size()];
      real W[tensor::W::size()];
    };
    struct ViscoElasticNeighborData {
      real w[tensor::w::size()];
    };
  }
}

#endif
