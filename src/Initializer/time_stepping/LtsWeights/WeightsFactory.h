#ifndef SEISSOL_LTSWEIGHTSFACTORY_H
#define SEISSOL_LTSWEIGHTSFACTORY_H

#include <stdexcept>
#include <sstream>
#include <memory>
#include "WeightsModels.h"


namespace seissol::initializers::time_stepping {

enum class NodeWeightModelTypes : int {
  ExponentialWeights = 0,
  ExponentialBalancedWeights,
  EncodedBalancedWeights,
  ExponentialBalancedWeightsWithBalancedMessaging,
  ExponentialBalancedWeightsWithMessageCount,
  Count
};

enum class EdgeWeightModelTypes : int {
  Naive = 0,
  ApproximateCommunication,
  ApproximateCommunicationWithBalancedMessaging,
  ApproximateCommunicationWithMessageCount,
  Count
};

inline void checkWeightModelCombinationIsAllowed(NodeWeightModelTypes nmod, EdgeWeightModelTypes emod) {
  if (nmod == NodeWeightModelTypes::ExponentialBalancedWeightsWithBalancedMessaging &&
      emod != EdgeWeightModelTypes::ApproximateCommunicationWithBalancedMessaging)
  {
    std::stringstream err;
    err << "Balanced Messaging strategies need to be paired";
    throw std::runtime_error(err.str());
  }

  if (nmod == NodeWeightModelTypes::ExponentialBalancedWeightsWithMessageCount && 
      emod != EdgeWeightModelTypes::ApproximateCommunicationWithMessageCount)
  {
    std::stringstream err;
    err << "Minimum Messaging strategies need to be paired";
    throw std::runtime_error(err.str());
  }

}

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

  checkWeightModelCombinationIsAllowed(nodeType, edgeType);

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
  case NodeWeightModelTypes::ExponentialBalancedWeightsWithBalancedMessaging: {
    nwm = new ExponentialBalancedWeightsWithBalancedMessaging(*lts.get());
    break;
  }
  case NodeWeightModelTypes::ExponentialBalancedWeightsWithMessageCount: {
    nwm = new ExponentialBalancedWeightsWithMessageCount(*lts.get());
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
  case EdgeWeightModelTypes::ApproximateCommunication: {
    ewm = new ApproximateCommunication(*lts.get());
    break;
  }
  case EdgeWeightModelTypes::ApproximateCommunicationWithBalancedMessaging: {
    ewm = new ApproximateCommunicationWithBalancedMessaging(*lts.get());
    break;
  }
  case EdgeWeightModelTypes::ApproximateCommunicationWithMessageCount:
  {
    ewm = new ApproximateCommunicationWithMessageCount(*lts.get());
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

#endif //SEISSOL_LTSWEIGHTSFACTORY_H
