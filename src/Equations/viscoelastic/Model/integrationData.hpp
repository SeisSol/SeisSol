#ifndef MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_
#define MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_

#include <generated_code/tensor.h>

#include <type_traits>
#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"

namespace seissol {
  namespace model {
    template<typename Config, int Mechanisms, std::enable_if_t<std::is_same_v<typename Config::MaterialT, ViscoElasticMaterial<Mechanisms>>, bool> = true>
    struct LocalSpecificData {
      typename Config::RealT E[tensor::E::size()];
      typename Config::RealT w[tensor::w::size()];
      typename Config::RealT W[tensor::W::size()];
    };

    template<typename Config, int Mechanisms, std::enable_if_t<std::is_same_v<typename Config::MaterialT, ViscoElasticMaterial<Mechanisms>>, bool> = true>
    struct NeighborSpecificData {
      typename Config::RealT w[tensor::w::size()];
    };
  }
}

#endif
