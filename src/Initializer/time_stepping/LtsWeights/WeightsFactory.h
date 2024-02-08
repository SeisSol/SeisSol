#ifndef SEISSOL_LTSWEIGHTSFACTORY_H
#define SEISSOL_LTSWEIGHTSFACTORY_H

#include <memory>
#include <sstream>
#include <stdexcept>

#include "WeightsModels.h"
#include "Initializer/Parameters/SeisSolParameters.h"

namespace seissol::initializer::time_stepping {

inline bool isLtsWeightsTypeAllowed(int id) {
  return ((id >= 0) && (id < static_cast<int>(parameters::LtsWeightsTypes::Count)));
}

inline parameters::LtsWeightsTypes convertLtsIdToType(int id) {
  if (isLtsWeightsTypeAllowed(id)) {
    return static_cast<parameters::LtsWeightsTypes>(id);
  }
  else {
    std::stringstream err;
    err << "provided LtsTWeightsType ID (" << id << ") is unknown";
    throw std::runtime_error(err.str());
  }
}

inline std::unique_ptr<LtsWeights> getLtsWeightsImplementation(parameters::LtsWeightsTypes type,
                                                        const LtsWeightsConfig& config,
                                                        seissol::SeisSol& seissolInstance) {
  switch (type) {
    case parameters::LtsWeightsTypes::ExponentialWeights : {
      return std::make_unique<ExponentialWeights>(config, seissolInstance);
    }
    case parameters::LtsWeightsTypes::ExponentialBalancedWeights : {
      return std::make_unique<ExponentialBalancedWeights>(config, seissolInstance);
    }
    case parameters::LtsWeightsTypes::EncodedBalancedWeights : {
      return std::make_unique<EncodedBalancedWeights>(config, seissolInstance);
    }
    default : {
      return std::unique_ptr<LtsWeights>(nullptr);
    }
  }
}

}

#endif //SEISSOL_LTSWEIGHTSFACTORY_H
