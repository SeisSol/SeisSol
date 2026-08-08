// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "Initializer/TimeStepping/ClusterSearch.h"

#include "Initializer/TimeStepping/ClusterHistogram.h"
#include "Initializer/TimeStepping/ClusterLadder.h"
#include "Initializer/TimeStepping/LtsWeights/LtsWeights.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <mpi.h>
#include <optional>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

namespace {
using time_stepping::computeGlobalCostOfClustering;
using time_stepping::computeMaxClusterIdAfterAutoMerge;
using time_stepping::enforceMaxClusterId;
} // namespace

ClusteringEvaluator::ClusteringEvaluator(const geometry::PumlMesh& mesh,
                                         parameters::BoundaryFormat boundaryFormat,
                                         const GlobalTimestep& timesteps,
                                         const std::vector<int>& cellCosts,
                                         std::vector<std::uint64_t> rate,
                                         bool smoothDuringSearch)
    : smoother_(mesh, boundaryFormat), timesteps_(&timesteps), cellCosts_(&cellCosts),
      rate_(std::move(rate)), smoothDuringSearch_(smoothDuringSearch),
      cellCount_(mesh.cells().size()), clusterIds_(mesh.cells().size(), 0) {}

std::vector<int> ClusteringEvaluator::binCells(double wiggleFactor) const {
  std::vector<int> clusterIds(cellCount_, 0);

  // Build the ladder once instead of walking the rate vector per cell.
  const auto ladder = ClusterLadder::forBinning(
      rate_, timesteps_->globalMinTimeStep, wiggleFactor, timesteps_->globalMaxTimeStep);

#pragma omp parallel for
  for (std::size_t cell = 0; cell < cellCount_; ++cell) {
    clusterIds[cell] = static_cast<int>(ladder.clusterOf(timesteps_->cellTimeStepWidths[cell]));
  }
  return clusterIds;
}

int ClusteringEvaluator::realize(double wiggleFactor) {
  int numberOfReductions = 0;
  auto lb = cache_.lower_bound(wiggleFactor);

  if (lb != cache_.end() && !(cache_.key_comp()(wiggleFactor, lb->first))) {
    clusterIds_ = lb->second;
  } else {
    // re-use best computed maxdiff enforcement available
    // reason that works: cf. Lukas' proof for cluster merging not violating maximum difference
    // we may generalize due to the fact that min(a, min(b,c)) = min(min(a,b), c) = min(min(a,c),
    // b), essentially establishing a partial ordering of clusterings, where A >= B iff
    // cluster(A[i]) >= cluster(B[i]) for all cells i. Thus: walking through the wiggle factors from
    // lower to higher will save a lot of reductions

    int cellchanges = 0;
    if (lb != cache_.end()) {
      // use the cache
      const auto newClusterIds = binCells(wiggleFactor);

#pragma omp parallel for reduction(+ : cellchanges)
      for (std::size_t cell = 0; cell < cellCount_; ++cell) {
        if (lb->second[cell] > newClusterIds[cell]) {
          ++cellchanges;
        }
        clusterIds_[cell] = std::min(lb->second[cell], newClusterIds[cell]);
      }
    } else {
      clusterIds_ = binCells(wiggleFactor);
      cellchanges = static_cast<int>(cellCount_);
    }
    if (smoothDuringSearch_) {
      MPI_Allreduce(MPI_IN_PLACE, &cellchanges, 1, MPI_INT, MPI_SUM, seissol::Mpi::mpi.comm());
      if (cellchanges > 0) {
        numberOfReductions = smoothCurrent();
      }
    }
    cache_[wiggleFactor] = clusterIds_;
  }

  return numberOfReductions;
}

int ClusteringEvaluator::smoothCurrent() {
  return smoother_.relax(clusterIds_, smoothingRule_, seissol::Mpi::mpi.comm());
}

int ClusteringEvaluator::globalMaxClusterId() const {
  int maxClusterId =
      clusterIds_.empty() ? 0 : *std::max_element(clusterIds_.begin(), clusterIds_.end());
  MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, Mpi::mpi.comm());
  return maxClusterId;
}

double ClusteringEvaluator::globalCost(double wiggleFactor) const {
  return computeGlobalCostOfClustering(clusterIds_,
                                       *cellCosts_,
                                       rate_,
                                       wiggleFactor,
                                       timesteps_->globalMinTimeStep,
                                       Mpi::mpi.comm());
}

ClusterHistogram ClusteringEvaluator::globalHistogram() const {
  auto histogram = ClusterHistogram::fromClustering(
      clusterIds_, *cellCosts_, static_cast<std::size_t>(globalMaxClusterId()) + 1);
  histogram.reduce(Mpi::mpi.comm());
  return histogram;
}

SearchResult GridLadderSearch::run(ClusteringEvaluator& evaluator,
                                   const SearchConstraints& constraints) {
  reductions_ = 0;

  if (constraints.autoMergeBaseline == parameters::AutoMergeCostBaseline::BestWiggleFactor) {
    // First compute wiggle factor without merging as baseline cost
    logInfo() << "Using best wiggle factor as baseline cost for auto merging.";
    logInfo() << "1. Compute best wiggle factor without merging clusters";
    const auto baseline = sweep(evaluator, constraints, std::nullopt, false);
    // Compute wiggle factor a second time with merging and using the previous cost as baseline
    logInfo() << "2. Compute best wiggle factor with merging clusters, using the previous cost "
                 "estimate as baseline";
    return sweep(evaluator, constraints, baseline.cost, constraints.autoMerge);
  }

  assert(constraints.autoMergeBaseline == parameters::AutoMergeCostBaseline::MaxWiggleFactor);
  return sweep(evaluator, constraints, std::nullopt, constraints.autoMerge);
}

SearchResult GridLadderSearch::sweep(ClusteringEvaluator& evaluator,
                                     const SearchConstraints& constraints,
                                     std::optional<double> baselineCost,
                                     bool autoMerge) {
  // Maps that keep track of number of clusters vs cost
  auto mapMaxClusterIdToLowestCost = std::map<int, double>{};
  auto maxMapClusterIdToBestWiggleFactor = std::map<int, double>{};

  const double minWiggleFactor = constraints.minWiggleFactor;
  const double maxWiggleFactor = 1.0;

  const double stepSizeWiggleFactor = constraints.wiggleFactorStepsize;
  const int numberOfStepsWiggleFactor =
      std::ceil((maxWiggleFactor - minWiggleFactor) / stepSizeWiggleFactor) + 1;

  auto computeWiggleFactor = [minWiggleFactor, stepSizeWiggleFactor, maxWiggleFactor](auto ith) {
    return std::min(minWiggleFactor + ith * stepSizeWiggleFactor, maxWiggleFactor);
  };

  auto totalWiggleFactorReductions = 0U;

  if (baselineCost) {
    logInfo() << "Baseline cost before cluster merging is" << *baselineCost;
  } else {
    // Compute baselineCost cost before wiggle factor and merging of clusters
    totalWiggleFactorReductions += evaluator.realize(maxWiggleFactor);
    baselineCost = evaluator.globalCost(maxWiggleFactor);
    logInfo() << "Baseline cost, without wiggle factor and cluster merging is" << *baselineCost;
  }
  assert(baselineCost);

  const double maxAdmissibleCost = constraints.allowedPerformanceLossRatio * *baselineCost;

  if (autoMerge) {
    logInfo() << "Maximal admissible cost after cluster merging is" << maxAdmissibleCost;
  }

  for (int i = 0; i < numberOfStepsWiggleFactor; ++i) {
    const double curWiggleFactor = computeWiggleFactor(i);
    totalWiggleFactorReductions += evaluator.realize(curWiggleFactor);

    // Note: Merging clusters does not invalidate invariance generated by enforceMaximumDifference()
    // This can be shown by enumerating all possible cases
    auto maxClusterIdToEnforce = constraints.maxClusterId;
    if (autoMerge) {
      const auto maxClusterIdAfterMerging =
          computeMaxClusterIdAfterAutoMerge(evaluator.clusterIds(),
                                            evaluator.cellCosts(),
                                            evaluator.rate(),
                                            maxAdmissibleCost,
                                            curWiggleFactor,
                                            evaluator.timesteps().globalMinTimeStep);
      maxClusterIdToEnforce = std::min(maxClusterIdAfterMerging, maxClusterIdToEnforce);
    }

    evaluator.mutableClusterIds() =
        enforceMaxClusterId(evaluator.clusterIds(), maxClusterIdToEnforce);
    const auto maxClusterId = evaluator.globalMaxClusterId();

    // Compute cost
    const double cost = evaluator.globalCost(curWiggleFactor);

    if (auto it = mapMaxClusterIdToLowestCost.find(maxClusterId);
        it == mapMaxClusterIdToLowestCost.end() || cost <= it->second) {
      maxMapClusterIdToBestWiggleFactor[maxClusterId] = curWiggleFactor;
      mapMaxClusterIdToLowestCost[maxClusterId] = cost;
    }
  }

  // Find best wiggle factor after merging of clusters
  // We compare against cost of baselineCost.
  int minAdmissibleMaxClusterId = std::numeric_limits<int>::max();
  if (autoMerge) {
    // When merging clusters, we want to find the minimum number of clusters with admissible
    // performance.
    bool foundAdmissibleMerge = false;
    for (const auto& [noOfClusters, cost] : mapMaxClusterIdToLowestCost) {
      if (cost <= maxAdmissibleCost) {
        foundAdmissibleMerge = true;
        minAdmissibleMaxClusterId = std::min(minAdmissibleMaxClusterId, noOfClusters);
        logDebug() << "Admissible. cluster:" << noOfClusters << ",cost" << cost
                   << "with wiggle factor" << maxMapClusterIdToBestWiggleFactor[noOfClusters];
      } else {
        logDebug() << "Not admissible. cluster:" << noOfClusters << ",cost" << cost
                   << "with wiggle factor" << maxMapClusterIdToBestWiggleFactor[noOfClusters];
      }
    }
    if (!foundAdmissibleMerge) {
      logError() << "Found no admissible wiggle factor with cluster merge. Aborting.";
    }
  } else {
    // Otherwise choose one with the smallest cost
    minAdmissibleMaxClusterId =
        std::min_element(mapMaxClusterIdToLowestCost.begin(),
                         mapMaxClusterIdToLowestCost.end(),
                         [](const auto& a, const auto& b) { return a.second < b.second; })
            ->first;
  }

  logInfo() << "Enforcing maximum difference when finding best wiggle factor took"
            << totalWiggleFactorReductions << "reductions.";
  reductions_ += static_cast<int>(totalWiggleFactorReductions);

  const auto bestWiggleFactor = maxMapClusterIdToBestWiggleFactor[minAdmissibleMaxClusterId];
  const auto bestCostEstimate = mapMaxClusterIdToLowestCost[minAdmissibleMaxClusterId];
  logInfo() << "The best wiggle factor is" << bestWiggleFactor << "with cost" << bestCostEstimate
            << "and" << minAdmissibleMaxClusterId + 1 << "time clusters";

  if (baselineCost > bestCostEstimate) {
    logInfo() << "Cost decreased" << (*baselineCost - bestCostEstimate) / *baselineCost * 100
              << "% with absolute cost decrease of" << *baselineCost - bestCostEstimate
              << "compared to the baseline";
  } else {
    logInfo() << "Cost increased" << (bestCostEstimate - *baselineCost) / *baselineCost * 100
              << "% with absolute cost increase of" << bestCostEstimate - *baselineCost
              << "compared to the baseline";
    logInfo() << "Note: Cost increased due to cluster merging!";
  }

  return SearchResult{bestWiggleFactor, minAdmissibleMaxClusterId, bestCostEstimate};
}

} // namespace seissol::initializer
