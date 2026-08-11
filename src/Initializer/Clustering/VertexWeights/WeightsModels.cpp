// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WeightsModels.h"

#include "GeneratedCode/init.h"
#include "Initializer/Clustering/ClusteringCost.h"
#include "Initializer/Clustering/VertexWeights/VertexWeightModel.h"

#include <cassert>
#include <cstddef>

namespace seissol::initializer {

void ExponentialWeights::setVertexWeights() {
  assert(ncon_ == 1 && "single constraint partitioning");
  const auto& clustering = this->clustering();

  for (std::size_t cell = 0; cell < clustering.cellCosts.size(); ++cell) {
    const auto factor = ratepow(clustering.ratios, clustering.clusterIds[cell], maxCluster());
    vertexWeights_[ncon_ * cell] = factor * clustering.cellCosts[cell];
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
  const auto& clustering = this->clustering();

  for (std::size_t cell = 0; cell < clustering.cellCosts.size(); ++cell) {
    const auto factor = ratepow(clustering.ratios, clustering.clusterIds[cell], maxCluster());
    vertexWeights_[ncon_ * cell] = factor * clustering.cellCosts[cell];

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
  // One constraint per cluster of the ladder that was actually chosen. Deriving it from the
  // configured rate vector, as this used to, is wrong as soon as a search picks a different
  // ladder.
  return static_cast<int>(clustering().clusterCount());
}

void EncodedBalancedWeights::setVertexWeights() {
  const auto& clustering = this->clustering();
  for (std::size_t cell = 0; cell < clustering.cellCosts.size(); ++cell) {
    for (int i = 0; i < ncon_; ++i) {
      vertexWeights_[ncon_ * cell + i] = 0;
    }
    vertexWeights_[ncon_ * cell + clustering.clusterIds[cell]] = clustering.cellCosts[cell];
  }
}

void EncodedBalancedWeights::setAllowedImbalances() {
  imbalances_.resize(ncon_);

  constexpr double MediumLtsWeightImbalance{1.05};
  for (int i = 0; i < ncon_; ++i) {
    imbalances_[i] = MediumLtsWeightImbalance;
  }
}
} // namespace seissol::initializer
