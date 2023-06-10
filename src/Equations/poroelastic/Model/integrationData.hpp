#ifndef MODEL_POROELASTIC_INTEGRATIONDATA_H_
#define MODEL_POROELASTIC_INTEGRATIONDATA_H_

#include <generated_code/tensor.h>
#include <type_traits>
#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"

namespace seissol {
  namespace model {

    template<typename Config, std::enable_if<std::is_same_v<Config::MaterialT, PoroElasticMaterial>, bool> = true>
    struct LocalSpecificData {
      Config::RealT sourceMatrix[seissol::tensor::ET::size()];
      Config::RealT G[Config::MaterialT::NumberOfQuantities];
      Config::RealT typicalTimeStepWidth;
      Config::RealT Zinv[Config::MaterialT::NumberOfQuantities][Config::ConvergenceOrder*Config::ConvergenceOrder];
    };

    template<typename Config, std::enable_if<std::is_same_v<Config::MaterialT, PoroElasticMaterial>, bool> = true>
    struct NeighborSpecificData {};
  }
}

#endif
