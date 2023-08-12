#ifndef MODEL_POROELASTIC_INTEGRATIONDATA_H_
#define MODEL_POROELASTIC_INTEGRATIONDATA_H_

#include <generated_code/tensor.h>
#include <type_traits>
#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"

namespace seissol {
  namespace model {

    template<typename Config>
    struct LocalSpecificData<Config, std::enable_if_t<std::is_same_v<typename Config::MaterialT, PoroElasticMaterial>>> {
      typename Config::RealT sourceMatrix[seissol::tensor::ET::size()];
      typename Config::RealT G[Config::MaterialT::NumberOfQuantities];
      typename Config::RealT typicalTimeStepWidth;
      typename Config::RealT Zinv[Config::MaterialT::NumberOfQuantities][Config::ConvergenceOrder*Config::ConvergenceOrder];
    };

    template<typename Config>
    struct NeighborSpecificData<Config, std::enable_if_t<std::is_same_v<typename Config::MaterialT, PoroElasticMaterial>>> {};
  }
}

#endif
