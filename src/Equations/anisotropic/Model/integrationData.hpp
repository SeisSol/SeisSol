#ifndef MODEL_ANISOTROPIC_INTEGRATIONDATA_H_
#define MODEL_ANISOTROPIC_INTEGRATIONDATA_H_

#include <type_traits>
#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"

namespace seissol {
  namespace model {

    template<typename Config, std::enable_if<std::is_same_v<Config::MaterialT, AnisotropicMaterial>, bool> = true>
    struct LocalSpecificData {
    };

    template<typename Config, std::enable_if<std::is_same_v<Config::MaterialT, AnisotropicMaterial>, bool> = true>
    struct NeighborSpecificData {
    };
  }
}

#endif
