// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_WEIGHTSFACTORY_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_WEIGHTSFACTORY_H_

#include <memory>
#include <sstream>
#include <stdexcept>

#include "Initializer/Parameters/SeisSolParameters.h"
#include "WeightsModels.h"

namespace seissol::initializer::time_stepping {

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

inline std::unique_ptr<LtsWeights> getLtsWeightsImplementation(parameters::LtsWeightsTypes type,
                                                               const LtsWeightsConfig& config,
                                                               seissol::SeisSol& seissolInstance) {
  switch (type) {
  case parameters::LtsWeightsTypes::ExponentialWeights: {
    return std::make_unique<ExponentialWeights>(config, seissolInstance);
  }
  case parameters::LtsWeightsTypes::ExponentialBalancedWeights: {
    return std::make_unique<ExponentialBalancedWeights>(config, seissolInstance);
  }
  case parameters::LtsWeightsTypes::EncodedBalancedWeights: {
    return std::make_unique<EncodedBalancedWeights>(config, seissolInstance);
  }
  default: {
    return std::unique_ptr<LtsWeights>(nullptr);
  }
  }
}

} // namespace seissol::initializer::time_stepping

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_WEIGHTSFACTORY_H_
