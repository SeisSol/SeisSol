#ifndef MODEL_ELASTIC_INTEGRATIONDATA_H_
#define MODEL_ELASTIC_INTEGRATIONDATA_H_

#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"
#include <type_traits>

namespace seissol {
  namespace model {

    template<typename Config, std::enable_if_t<std::is_same_v<typename Config::MaterialT, ElasticMaterial>, bool> = true>
    struct LocalSpecificData {
    };

    template<typename Config, std::enable_if_t<std::is_same_v<typename Config::MaterialT, ElasticMaterial>, bool> = true>
    struct NeighborSpecificData {
    };
  }
}

#endif
