// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_VERTEXWEIGHTMODEL_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_VERTEXWEIGHTMODEL_H_

#include "Initializer/Clustering/Clustering.h"

#include <cstddef>
#include <vector>

namespace seissol::initializer {

/// Turns a clustering into the vertex weights the graph partitioner is fed with.
///
/// This is the only part of the LTS setup that is polymorphic, and the only part that knows
/// about ParMETIS conventions. It reads a ClusteringResult and produces weights; it has no say
/// in how the clustering was arrived at.
class VertexWeightModel {
  public:
  VertexWeightModel() = default;
  virtual ~VertexWeightModel() = default;
  VertexWeightModel(const VertexWeightModel&) = delete;
  VertexWeightModel& operator=(const VertexWeightModel&) = delete;
  VertexWeightModel(VertexWeightModel&&) = delete;
  VertexWeightModel& operator=(VertexWeightModel&&) = delete;

  /// `clustering` must outlive the call, not the object.
  void build(const ClusteringResult& clustering);

  [[nodiscard]] const int* vertexWeights() const;
  [[nodiscard]] const double* imbalances() const;
  [[nodiscard]] int nWeightsPerVertex() const;

  protected:
  virtual int evaluateNumberOfConstraints() = 0;
  virtual void setVertexWeights() = 0;
  virtual void setAllowedImbalances() = 0;

  [[nodiscard]] const ClusteringResult& clustering() const { return *clustering_; }
  /// Update factor exponent of the coarsest cluster; weights are scaled relative to it so
  /// that they stay integral.
  [[nodiscard]] std::size_t maxCluster() const { return clustering_->ratios.size(); }

  const ClusteringResult* clustering_{nullptr};
  std::vector<int> vertexWeights_;
  std::vector<double> imbalances_;
  int ncon_{0};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_VERTEXWEIGHTMODEL_H_
