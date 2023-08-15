#ifndef MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_
#define MODEL_VISCOELASTIC2_INTEGRATIONDATA_H_

#include <generated_code/tensor.h>

#include "Common/cellconfigconv.hpp"
#include <type_traits>
#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"
#include "Common/configtensor.hpp"

namespace seissol {
  namespace model {
    template<typename Config>
    struct LocalSpecificData<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevskiAnelastic>> {
      typename Config::RealT E[ConfigConstants<Config>::TensorSizeE];
      typename Config::RealT w[ConfigConstants<Config>::TensorSizew];
      typename Config::RealT W[ConfigConstants<Config>::TensorSizeW];
    };

    template<typename Config>
    struct NeighborSpecificData<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevskiAnelastic>> {
      typename Config::RealT w[ConfigConstants<Config>::TensorSizew];
    };
  }
}

#endif
