#ifndef SEISSOL_LTSWEIGHTSFACTORY_H
#define SEISSOL_LTSWEIGHTSFACTORY_H

#include <stdexcept>
#include <sstream>
#include <memory>
#include "WeightsModels.h"


namespace seissol::initializers::time_stepping {

enum class LtsWeightsTypes: int {
  ExponentialWeights = 0,
  ExponentialBalancedWeights,
  EncodedBalancedWeights,
  Count
};

inline bool isLtsWeightsTypeAllowed(int id) {
  return ((id >= 0) && (id < static_cast<int>(LtsWeightsTypes::Count)));
}

inline LtsWeightsTypes convertLtsIdToType(int id) {
  if (isLtsWeightsTypeAllowed(id)) {
    return static_cast<LtsWeightsTypes>(id);
  }
  else {
    std::stringstream err;
    err << "provided LtsTWeightsType ID (" << id << ") is unknown";
    throw std::runtime_error(err.str());
  }
}

inline std::unique_ptr<LtsWeights> getLtsWeightsImplementation(LtsWeightsTypes type,
                                                        const LtsWeightsConfig& config,
                                                        const LtsParameters* ltsParameters) {
  switch (type) {
    case LtsWeightsTypes::ExponentialWeights : {
      return std::make_unique<ExponentialWeights>(config, ltsParameters);
    }
    case LtsWeightsTypes::ExponentialBalancedWeights : {
      return std::make_unique<ExponentialBalancedWeights>(config, ltsParameters);
    }
    case LtsWeightsTypes::EncodedBalancedWeights : {
      return std::make_unique<EncodedBalancedWeights>(config, ltsParameters);
    }
    default : {
      return std::unique_ptr<LtsWeights>(nullptr);
    }
  }
}

}

#endif //SEISSOL_LTSWEIGHTSFACTORY_H
