// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_LADDERSEARCH_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_LADDERSEARCH_H_

#include "Initializer/Clustering/ClusterCostModel.h"
#include "Initializer/Clustering/ClusteringEvaluator.h"
#include "Initializer/Parameters/LtsParameters.h"

#include <cstdint>
#include <vector>

namespace seissol::initializer {

struct SearchConstraints {
  double minWiggleFactor{1.0};
  double wiggleFactorStepsize{0.01};
  int maxClusterId{0};
  bool autoMerge{false};
  double allowedPerformanceLossRatio{1.0};
  parameters::AutoMergeCostBaseline autoMergeBaseline{
      parameters::AutoMergeCostBaseline::MaxWiggleFactor};
  /// Only the lattice search reads these: the legacy search does not choose the ladder, so
  /// it has no way to trade cluster count against fill.
  ClusterCostModel costModel{};
  std::uint64_t maxRatio{0};
};

struct SearchResult {
  double wiggleFactor{1.0};
  int maxClusterId{0};
  double cost{0.0};
  /// The ladder the search settled on, **complete**: one ratio per cluster boundary, with no
  /// repetition of the last entry and no terminator.
  ///
  /// The abbreviated form that `ClusteredLTS` accepts exists only between the parameter file
  /// and the search. Every search hands back a complete ladder, so that nothing downstream has
  /// to guess which of the two conventions a vector follows.
  std::vector<std::uint64_t> ratios;
};

/// Realizes candidate clusterings and scores them.
///
/// This owns everything mesh-derived -- the smoother, the cell timesteps, the cell costs --
/// so that a search implementation is a pure function of this object plus the constraints,
/// and needs no access to LtsWeights or the partitioner-facing weight models.
///
/// Realized clusterings are cached by wiggle factor. Sweeping wiggle factors upwards lets
/// each candidate start its smoothing fixed point from the previous, already-smoothed
/// clustering, because the cluster ids are monotone in the wiggle factor. That trick is
/// what keeps the number of reductions manageable, so the cache is part of the contract
/// rather than an optimization detail.
/// Picks a ladder. Implementations differ in how they explore, not in what they optimize.
class LadderSearch {
  public:
  LadderSearch() = default;
  virtual ~LadderSearch() = default;
  LadderSearch(const LadderSearch&) = delete;
  LadderSearch& operator=(const LadderSearch&) = delete;
  LadderSearch(LadderSearch&&) = delete;
  LadderSearch& operator=(LadderSearch&&) = delete;

  virtual SearchResult run(ClusteringEvaluator& evaluator,
                           const SearchConstraints& constraints) = 0;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_LADDERSEARCH_H_
