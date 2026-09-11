// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERING_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERING_H_

#include "Geometry/PUMLReader.h"
#include "Initializer/FaceMap.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer {

struct ClusteringConfig {
  parameters::BoundaryFormat boundaryFormat;
  std::vector<std::uint64_t> rate;
  std::uint64_t vertexWeightElement{};
  std::uint64_t vertexWeightDynamicRupture{};
  std::uint64_t vertexWeightFreeSurfaceWithGravity{};
  const FaceMap* faceMap{nullptr};
};

/// Everything the clustering decided.
///
/// This is the whole interface between deciding the clustering and using it: the partitioner
/// weights are derived from it, the mesh reader carries it, and the ladder plus the wiggle
/// factor travel through it to ClusterLayout. Nothing downstream may re-derive any of it from
/// the parameter file, because a search is allowed to have changed both.
struct ClusteringResult {
  /// Cluster id per cell, in the order the mesh had before partitioning.
  std::vector<std::size_t> clusterIds;
  /// The ladder, complete: one ratio per cluster boundary, never abbreviated.
  std::vector<std::uint64_t> ratios;
  double wiggleFactor{1.0};
  /// Per-cell cost the clustering was scored with; also the base of the vertex weights.
  std::vector<std::uint64_t> cellCosts;
  GlobalTimestep timesteps;

  [[nodiscard]] std::size_t clusterCount() const { return ratios.size() + 1; }
};

/// Decides which cell belongs to which time cluster.
///
/// Everything about *how* that decision is reached -- the ladder, the histograms, the
/// smoothing, the search -- lives behind this class. What comes out is a ClusteringResult;
/// turning that into partitioner weights is a separate concern, see VertexWeightModel.
class Clustering {
  public:
  Clustering(const ClusteringConfig& config, seissol::SeisSol& seissolInstance);

  const ClusteringResult& compute(const geometry::PumlMesh& meshTopology,
                                  const geometry::PumlMesh& meshGeometry);

  [[nodiscard]] const ClusteringResult& result() const { return result_; }

  private:
  [[nodiscard]] std::vector<std::uint64_t>
      computeCostsPerTimestep(const geometry::PumlMesh& mesh) const;

  seissol::SeisSol& seissolInstance_;
  std::vector<std::uint64_t> rate_;
  std::uint64_t vertexWeightElement_{};
  std::uint64_t vertexWeightDynamicRupture_{};
  std::uint64_t vertexWeightFreeSurfaceWithGravity_{};
  parameters::BoundaryFormat boundaryFormat_;
  const FaceMap* faceMap_{nullptr};

  ClusteringResult result_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERING_H_
