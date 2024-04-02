#ifndef MODEL_POROELASTIC_INTEGRATIONDATA_H_
#define MODEL_POROELASTIC_INTEGRATIONDATA_H_

#include "datastructures.hpp"
#include <Common/constants.hpp>
#include <generated_code/tensor.h>

namespace seissol {
  namespace model {

    struct PoroelasticLocalData {
      real sourceMatrix[seissol::tensor::ET::size()];
      real G[seissol::model::Material_t::NumberOfQuantities];
      real typicalTimeStepWidth;
      real Zinv[seissol::model::Material_t::NumberOfQuantities][ConvergenceOrder*ConvergenceOrder];
    };
    struct PoroelasticNeighborData {};
  }
}

#endif
