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
#include <cstdint>
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

  /// Weights as the graph partitioner wants them: one block of `nWeightsPerVertex()` per cell.
  /// Unsigned and 64 bit wide because a weight is `cellCost * updateFactor`, which overflows a
  /// 32-bit int for deep ladders; PUML widens whatever it is given to `unsigned long` anyway.
  [[nodiscard]] const std::uint64_t* vertexWeights() const;
  [[nodiscard]] const double* imbalances() const;
  /// Kept signed because PUML's `setVertexWeights` takes an `int` count; everything on this
  /// side of that call uses `ncon_`.
  [[nodiscard]] int nWeightsPerVertex() const;

  protected:
  virtual std::size_t evaluateNumberOfConstraints() = 0;
  virtual void setVertexWeights() = 0;
  virtual void setAllowedImbalances() = 0;

  [[nodiscard]] const ClusteringResult& clustering() const { return *clustering_; }
  /// Update factor exponent of the coarsest cluster; weights are scaled relative to it so
  /// that they stay integral.
  [[nodiscard]] std::size_t maxCluster() const { return clustering_->ratios.size(); }

  const ClusteringResult* clustering_{nullptr};
  std::vector<std::uint64_t> vertexWeights_;
  std::vector<double> imbalances_;
  std::size_t ncon_{0};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_VERTEXWEIGHTMODEL_H_
