// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERINGEVALUATOR_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERINGEVALUATOR_H_

#include "Geometry/PUMLReader.h"
#include "Initializer/Clustering/ClusterHistogram.h"
#include "Initializer/Clustering/ClusterLadder.h"
#include "Initializer/Clustering/ClusterSmoother.h"
#include "Initializer/Clustering/TimestepHistogram.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"

#include <cstdint>
#include <functional>
#include <map>
#include <vector>

namespace seissol::initializer {

class ClusteringEvaluator {
  public:
  /// `timesteps` and `cellCosts` are referenced, not copied: at the meshes this runs on,
  /// duplicating a per-cell array costs more setup memory than the class is worth. Both must
  /// outlive the evaluator.
  ClusteringEvaluator(const geometry::PumlMesh& mesh,
                      parameters::BoundaryFormat boundaryFormat,
                      const FaceMap& faceMap,
                      const GlobalTimestep& timesteps,
                      const std::vector<std::uint64_t>& cellCosts,
                      std::vector<std::uint64_t> rate,
                      bool smoothDuringSearch);

  /// Bins the cells for `wiggleFactor` and smooths the result. Returns the number of
  /// demotions performed, which is zero when the candidate was served from the cache.
  int realize(double wiggleFactor);

  /// Same, for an explicitly given *complete* ladder -- one ratio per cluster boundary, with
  /// no repetition of the last entry. Falls back to the cached path when the ladder happens to
  /// be the configured one; other ladders are realized from scratch, because the monotonicity
  /// argument behind the cache holds along the wiggle factor at a fixed ladder, not across
  /// ladders.
  int realize(const std::vector<std::uint64_t>& ratios, double wiggleFactor);

  /// The configured rate vector expanded to a complete ladder at `wiggleFactor`.
  [[nodiscard]] std::vector<std::uint64_t> configuredRatios(double wiggleFactor) const;

  /// Smooths the current clustering without touching the cache.
  int smoothCurrent();

  [[nodiscard]] const std::vector<std::size_t>& clusterIds() const { return clusterIds_; }
  std::vector<std::size_t>& mutableClusterIds() { return clusterIds_; }

  /// Largest cluster id present on any rank.
  [[nodiscard]] std::size_t globalMaxClusterId() const;

  /// Exact per-cell cost of the current clustering, summed over all ranks.
  [[nodiscard]] double globalCost(double wiggleFactor) const;
  [[nodiscard]] double globalCost(const std::vector<std::uint64_t>& ratios,
                                  double wiggleFactor) const;

  /// Cell weight over `floor(cellTimestep / (wiggleFactor * minimumTimestep))`, reduced over
  /// all ranks. One collective; afterwards any ladder on that base timestep can be scored
  /// without touching the mesh again.
  [[nodiscard]] TimestepHistogram timestepHistogram(double wiggleFactor,
                                                    std::size_t maxIndex) const;

  /// Per-cluster weights of the current clustering, summed over all ranks. The substrate
  /// for cost models that do not want to walk the cells again.
  [[nodiscard]] ClusterHistogram globalHistogram() const;

  [[nodiscard]] const GlobalTimestep& timesteps() const { return *timesteps_; }
  [[nodiscard]] const std::vector<std::uint64_t>& cellCosts() const { return *cellCosts_; }
  [[nodiscard]] const std::vector<std::uint64_t>& rate() const { return rate_; }

  private:
  [[nodiscard]] ClusterLadder configuredLadder(double wiggleFactor) const;
  std::vector<std::size_t> binCells(const ClusterLadder& ladder) const;

  ClusterSmoother smoother_;
  SmoothingRule smoothingRule_{};
  const GlobalTimestep* timesteps_;
  const std::vector<std::uint64_t>* cellCosts_;
  std::vector<std::uint64_t> rate_;
  bool smoothDuringSearch_;

  std::size_t cellCount_{0};
  std::vector<std::size_t> clusterIds_;
  // maps wiggle factor to the realized, uncapped clustering
  std::map<double, std::vector<std::size_t>, std::greater<>> cache_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERINGEVALUATOR_H_
