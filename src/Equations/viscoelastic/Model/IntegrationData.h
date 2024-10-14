#ifndef VISCOELASTIC_INTEGRATIONDATA_H_
#define VISCOELASTIC_INTEGRATIONDATA_H_

#include "generated-code/tensor.h"

namespace seissol::model {

struct ViscoElasticLocalData {
  real sourceMatrix[seissol::tensor::ET::size()];
};
struct ViscoElasticNeighborData {};
} // namespace seissol::model

#endif
