#ifndef MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_
#define MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_

#include <generated_code/tensor.h>

#include <type_traits>
#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"

namespace seissol {
  namespace model {
    template<typename Config>
    struct LocalSpecificData<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevskiAnelastic>> {
      typename Config::RealT E[tensor::E::size()];
      typename Config::RealT w[tensor::w::size()];
      typename Config::RealT W[tensor::W::size()];
    };

    template<typename Config>
    struct NeighborSpecificData<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevskiAnelastic>> {
      typename Config::RealT w[tensor::w::size()];
    };
  }
}

#endif
