#ifndef MODEL_POROELASTIC_INTEGRATIONDATA_H_
#define MODEL_POROELASTIC_INTEGRATIONDATA_H_

#include "Common/constants.hpp"
#include "datastructures.hpp"
#include "generated_code/tensor.h"

namespace seissol::model {

struct PoroelasticLocalData {
  real sourceMatrix[seissol::tensor::ET::size()];
  real G[PoroElasticMaterial::NumberOfQuantities];
  real typicalTimeStepWidth;
  real Zinv[PoroElasticMaterial::NumberOfQuantities][ConvergenceOrder * ConvergenceOrder];
};
struct PoroelasticNeighborData {};

} // namespace seissol::model

#endif
