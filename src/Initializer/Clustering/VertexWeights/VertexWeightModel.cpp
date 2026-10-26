// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Clustering/VertexWeights/VertexWeightModel.h"

#include "Initializer/Clustering/Clustering.h"

#include <cassert>
#include <cstdint>

namespace seissol::initializer {

void VertexWeightModel::build(const ClusteringResult& clustering) {
  clustering_ = &clustering;

  ncon_ = evaluateNumberOfConstraints();
  assert(ncon_ > 0 && "the number of constraints has to be positive");

  vertexWeights_.assign(clustering.clusterIds.size() * ncon_, 0);

  setVertexWeights();
  setAllowedImbalances();
}

const std::uint64_t* VertexWeightModel::vertexWeights() const {
  assert(!vertexWeights_.empty() && "vertex weights are not initialized");
  return vertexWeights_.data();
}

const double* VertexWeightModel::imbalances() const {
  assert(!imbalances_.empty() && "weight imbalances are not initialized");
  return imbalances_.data();
}

int VertexWeightModel::nWeightsPerVertex() const {
  assert(ncon_ > 0 && "num. constraints has not been initialized yet");
  return static_cast<int>(ncon_);
}

} // namespace seissol::initializer
