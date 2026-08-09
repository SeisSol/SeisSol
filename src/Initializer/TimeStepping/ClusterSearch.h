// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERSEARCH_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERSEARCH_H_

#include "Geometry/PUMLReader.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/TimeStepping/ClusterCostModel.h"
#include "Initializer/TimeStepping/ClusterHistogram.h"
#include "Initializer/TimeStepping/ClusterSmoother.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "Initializer/TimeStepping/TimestepHistogram.h"

#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <vector>

namespace seissol::initializer {

/// What the search is allowed to do, straight from the parameter file.
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
  /// The ladder the search settled on, in the user-facing (possibly abbreviated) form. A
  /// search that only tunes the wiggle factor hands back what it was given; one that also
  /// chooses the ladder hands back something else, and that is what has to reach
  /// ClusterLayout.
  std::vector<std::uint64_t> rates;
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
class ClusteringEvaluator {
  public:
  /// `timesteps` and `cellCosts` are referenced, not copied: at the meshes this runs on,
  /// duplicating a per-cell array costs more setup memory than the class is worth. Both must
  /// outlive the evaluator.
  ClusteringEvaluator(const geometry::PumlMesh& mesh,
                      parameters::BoundaryFormat boundaryFormat,
                      const FaceMap& faceMap,
                      const GlobalTimestep& timesteps,
                      const std::vector<int>& cellCosts,
                      std::vector<std::uint64_t> rate,
                      bool smoothDuringSearch);

  /// Bins the cells for `wiggleFactor` and smooths the result. Returns the number of
  /// demotions performed, which is zero when the candidate was served from the cache.
  int realize(double wiggleFactor);

  /// Same, for a ladder other than the configured one. Only the configured rate vector is
  /// cached: the monotonicity argument behind the cache holds along the wiggle factor at a
  /// fixed ladder, not across ladders.
  int realize(const std::vector<std::uint64_t>& rates, double wiggleFactor);

  /// Smooths the current clustering without touching the cache.
  int smoothCurrent();

  [[nodiscard]] const std::vector<int>& clusterIds() const { return clusterIds_; }
  std::vector<int>& mutableClusterIds() { return clusterIds_; }

  /// Largest cluster id present on any rank.
  [[nodiscard]] int globalMaxClusterId() const;

  /// Exact per-cell cost of the current clustering, summed over all ranks.
  [[nodiscard]] double globalCost(double wiggleFactor) const;
  [[nodiscard]] double globalCost(const std::vector<std::uint64_t>& rates,
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
  [[nodiscard]] const std::vector<int>& cellCosts() const { return *cellCosts_; }
  [[nodiscard]] const std::vector<std::uint64_t>& rate() const { return rate_; }

  private:
  std::vector<int> binCells(const std::vector<std::uint64_t>& rates, double wiggleFactor) const;

  ClusterSmoother smoother_;
  SmoothingRule smoothingRule_{};
  const GlobalTimestep* timesteps_;
  const std::vector<int>* cellCosts_;
  std::vector<std::uint64_t> rate_;
  bool smoothDuringSearch_;

  std::size_t cellCount_{0};
  std::vector<int> clusterIds_;
  // maps wiggle factor to the realized, uncapped clustering
  std::map<double, std::vector<int>, std::greater<>> cache_;
};

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

/// The wiggle factor grid sweep with a top-merge scan per candidate, as SeisSol has always
/// done it. Kept as the default so that existing parameter files keep their behavior.
class GridLadderSearch : public LadderSearch {
  public:
  SearchResult run(ClusteringEvaluator& evaluator, const SearchConstraints& constraints) override;

  /// Total number of demotions performed while searching, for reporting.
  [[nodiscard]] int reductions() const { return reductions_; }

  private:
  SearchResult sweep(ClusteringEvaluator& evaluator,
                     const SearchConstraints& constraints,
                     std::optional<double> baselineCost,
                     bool autoMerge);

  int reductions_{0};
};

/// Chooses the ladder as well as the wiggle factor.
///
/// For every wiggle factor on the same grid the legacy search uses, the cell timesteps are
/// binned by update factor and the cheapest divisibility chain is found by dynamic
/// programming. Candidate evaluation therefore costs one histogram reduction plus local
/// work, instead of a smoothing fixed point per candidate; only the winner is realized.
///
/// With the default cost model this degenerates to "as many clusters as possible", because
/// pure update counting never charges for a cluster. Setting a launch cost or a fill
/// threshold is what makes coarser ladders competitive.
class LatticeDpSearch : public LadderSearch {
  public:
  SearchResult run(ClusteringEvaluator& evaluator, const SearchConstraints& constraints) override;

  [[nodiscard]] int reductions() const { return reductions_; }

  private:
  int reductions_{0};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_CLUSTERSEARCH_H_
