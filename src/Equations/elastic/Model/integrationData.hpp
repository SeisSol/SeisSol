#ifndef MODEL_ELASTIC_INTEGRATIONDATA_H_
#define MODEL_ELASTIC_INTEGRATIONDATA_H_

#include "Common/cellconfigconv.hpp"
#include "Model/common_datastructures.hpp"
#include "datastructures.hpp"
#include <type_traits>

namespace seissol {
  namespace model {

    template<typename Config>
    struct LocalSpecificData<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevski>> {
    };

    template<typename Config>
    struct NeighborSpecificData<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevski>> {
    };
  }
}

#endif
