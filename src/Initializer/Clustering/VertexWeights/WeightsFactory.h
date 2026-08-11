// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_WEIGHTSFACTORY_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_WEIGHTSFACTORY_H_

#include "Initializer/Clustering/VertexWeights/WeightsModels.h"
#include "Initializer/Parameters/SeisSolParameters.h"

#include <memory>
#include <sstream>
#include <stdexcept>

namespace seissol::initializer {

inline bool isLtsWeightsTypeAllowed(int id) {
  return ((id >= 0) && (id < static_cast<int>(parameters::LtsWeightsTypes::Count)));
}

inline parameters::LtsWeightsTypes convertLtsIdToType(int id) {
  if (isLtsWeightsTypeAllowed(id)) {
    return static_cast<parameters::LtsWeightsTypes>(id);
  } else {
    std::stringstream err;
    err << "provided LtsTWeightsType ID (" << id << ") is unknown";
    throw std::runtime_error(err.str());
  }
}

inline std::unique_ptr<VertexWeightModel> getVertexWeightModel(parameters::LtsWeightsTypes type) {
  switch (type) {
  case parameters::LtsWeightsTypes::ExponentialWeights: {
    return std::make_unique<ExponentialWeights>();
  }
  case parameters::LtsWeightsTypes::ExponentialBalancedWeights: {
    return std::make_unique<ExponentialBalancedWeights>();
  }
  case parameters::LtsWeightsTypes::EncodedBalancedWeights: {
    return std::make_unique<EncodedBalancedWeights>();
  }
  default: {
    return std::unique_ptr<VertexWeightModel>(nullptr);
  }
  }
}

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_WEIGHTSFACTORY_H_
