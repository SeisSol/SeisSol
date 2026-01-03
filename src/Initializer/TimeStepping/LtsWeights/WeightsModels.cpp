// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WeightsModels.h"

#include "GeneratedCode/init.h"
#include "Initializer/TimeStepping/LtsWeights/LtsWeights.h"

#include <cassert>
#include <cstddef>

namespace seissol::initializer::time_stepping {

void ExponentialWeights::setVertexWeights() {
  assert(ncon_ == 1 && "single constraint partitioning");
  const int maxCluster =
      getCluster(details_.globalMaxTimeStep, details_.globalMinTimeStep, wiggleFactor, rate_);

  for (std::size_t cell = 0; cell < cellCosts_.size(); ++cell) {
    const auto factor = ratepow(rate_, clusterIds_[cell], maxCluster);
    vertexWeights_[ncon_ * cell] = factor * cellCosts_[cell];
  }
}

void ExponentialWeights::setAllowedImbalances() {
  assert(ncon_ == 1 && "single constraint partitioning");
  imbalances_.resize(ncon_);

  constexpr double TinyLtsWeightImbalance{1.01};
  imbalances_[0] = TinyLtsWeightImbalance;
}

void ExponentialBalancedWeights::setVertexWeights() {
  assert(ncon_ == 2 && "binary constaints partitioning");
  const int maxCluster =
      getCluster(details_.globalMaxTimeStep, details_.globalMinTimeStep, wiggleFactor, rate_);

  for (std::size_t cell = 0; cell < cellCosts_.size(); ++cell) {
    const auto factor = ratepow(rate_, clusterIds_[cell], maxCluster);
    vertexWeights_[ncon_ * cell] = factor * cellCosts_[cell];

    constexpr int MemoryWeight{1};
    vertexWeights_[ncon_ * cell + 1] = MemoryWeight;
  }
}

void ExponentialBalancedWeights::setAllowedImbalances() {
  assert(ncon_ == 2 && "binary constaints partitioning");
  imbalances_.resize(ncon_);

  constexpr double TinyLtsWeightImbalance{1.01};
  imbalances_[0] = TinyLtsWeightImbalance;

  constexpr double MediumLtsMemoryImbalance{1.05};
  imbalances_[1] = MediumLtsMemoryImbalance;
}

int EncodedBalancedWeights::evaluateNumberOfConstraints() {
  const int maxCluster =
      getCluster(details_.globalMaxTimeStep, details_.globalMinTimeStep, wiggleFactor, rate_);
  return maxCluster + 1;
}

void EncodedBalancedWeights::setVertexWeights() {
  for (std::size_t cell = 0; cell < cellCosts_.size(); ++cell) {
    for (int i = 0; i < ncon_; ++i) {
      vertexWeights_[ncon_ * cell + i] = 0;
    }
    vertexWeights_[ncon_ * cell + clusterIds_[cell]] = cellCosts_[cell];
  }
}

void EncodedBalancedWeights::setAllowedImbalances() {
  imbalances_.resize(ncon_);

  constexpr double MediumLtsWeightImbalance{1.05};
  for (int i = 0; i < ncon_; ++i) {
    imbalances_[i] = MediumLtsWeightImbalance;
  }
}
} // namespace seissol::initializer::time_stepping
