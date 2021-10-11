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

bool isLtsWeightsTypeAllowed(int id) {
  return ((id >= 0) && (id < static_cast<int>(LtsWeightsTypes::Count)));
}

LtsWeightsTypes convertLtsIdToType(int id) {
  if (isLtsWeightsTypeAllowed(id)) {
    return static_cast<LtsWeightsTypes>(id);
  }
  else {
    std::stringstream err;
    err << "provided LtsTWeightsType ID (" << id << ") is unknown";
    throw std::runtime_error(err.str());
  }
}

std::unique_ptr<LtsWeights> getLtsWeightsImplementation(LtsWeightsTypes type, const LtsWeightsConfig& config) {
  switch (type) {
    case LtsWeightsTypes::ExponentialWeights : {
      return std::make_unique<ExponentialWeights>(config);
    }
    case LtsWeightsTypes::ExponentialBalancedWeights : {
      return std::make_unique<ExponentialBalancedWeights>(config);
    }
    case LtsWeightsTypes::EncodedBalancedWeights : {
      return std::make_unique<EncodedBalancedWeights>(config);
    }
    default : {
      return std::unique_ptr<LtsWeights>(nullptr);
    }
  }
}

}

#endif //SEISSOL_LTSWEIGHTSFACTORY_H
