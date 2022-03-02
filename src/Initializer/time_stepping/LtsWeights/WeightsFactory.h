#ifndef SEISSOL_LTSWEIGHTSFACTORY_H
#define SEISSOL_LTSWEIGHTSFACTORY_H

#include "LtsWeights.h"
#include "WeightsModels.h"
#include <memory>
#include <sstream>
#include <stdexcept>

namespace seissol::initializers::time_stepping {

enum class NodeWeightModelTypes : int {
  ExponentialWeights = 0,
  ExponentialBalancedWeights,
  EncodedBalancedWeights,
  Count
};

enum class EdgeWeightModelTypes : int {
  Naive = 0,
  PenalizeBetweenClusters,
  ApproximateCommunication,
  Count
};

bool isNodeWeightModelTypeAllowed(int id) {
  return ((id >= 0) && (id < static_cast<int>(NodeWeightModelTypes::Count)));
}

NodeWeightModelTypes convertNodeWeightModelTypeIdToType(int id) {
  if (isNodeWeightModelTypeAllowed(id)) {
    return static_cast<NodeWeightModelTypes>(id);
  } else {
    std::stringstream err;
    err << "provided NodeWeightModelType ID (" << id << ") is unknown";
    throw std::runtime_error(err.str());
  }
}

bool isEdgeWeightModelTypeAllowed(int id) {
  return ((id >= 0) && (id < static_cast<int>(EdgeWeightModelTypes::Count)));
}

EdgeWeightModelTypes convertEdgeWeightModelTypeIdToType(int id) {
  if (isEdgeWeightModelTypeAllowed(id)) {
    return static_cast<EdgeWeightModelTypes>(id);
  } else {
    std::stringstream err;
    err << "provided EdgeWeightModelType ID (" << id << ") is unknown";
    throw std::runtime_error(err.str());
  }
}

std::unique_ptr<LtsWeights> getLtsWeightsImplementation(NodeWeightModelTypes nodeType,
                                                        EdgeWeightModelTypes edgeType,
                                                        const LtsWeightsConfig& config) {
  std::unique_ptr<LtsWeights> lts = std::make_unique<LtsWeights>(config);
  NodeWeightModel* nwm = nullptr;
  EdgeWeightModel* ewm = nullptr;

  switch (nodeType) {
  case NodeWeightModelTypes::ExponentialWeights: {
    nwm = new ExponentialWeights(*lts.get());
    break;
  }
  case NodeWeightModelTypes::ExponentialBalancedWeights: {
    nwm = new ExponentialBalancedWeights(*lts.get());
    break;
  }
  case NodeWeightModelTypes::EncodedBalancedWeights: {
    nwm = new EncodedBalancedWeights(*lts.get());
    break;
  }
  default: {
    break;
  }
  }

  switch (edgeType) {
  case EdgeWeightModelTypes::Naive: {
    ewm = new Naive(*lts.get());
    break;
  }
  case EdgeWeightModelTypes::PenalizeBetweenClusters: {
    ewm = new PenalizeBetweenClusters(*lts.get());
    break;
  }
  case EdgeWeightModelTypes::ApproximateCommunication: {
    ewm = new ApproximateCommunication(*lts.get());
    break;
  }
  default: {
    break;
  }
  }

  lts->addWeightModels(nwm, ewm);
  return lts;
}

} // namespace seissol::initializers::time_stepping

#endif
