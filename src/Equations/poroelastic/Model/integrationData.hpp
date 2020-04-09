#ifndef MODEL_POROELASTIC_INTEGRATIONDATA_H_
#define MODEL_POROELASTIC_INTEGRATIONDATA_H_

namespace seissol {
  namespace model {

    struct PoroelasticLocalData {
      real sourceMatrix[seissol::tensor::ET::size()];
      real Zinv[NUMBER_OF_QUANTITIES*CONVERGENCE_ORDER*CONVERGENCE_ORDER];
    };
    struct PoroelasticNeighborData {};
  }
}

#endif
