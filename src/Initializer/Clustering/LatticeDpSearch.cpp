// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "Initializer/Clustering/LatticeDpSearch.h"

#include "Initializer/Clustering/ClusteringEvaluator.h"
#include "Initializer/Clustering/LadderOptimizer.h"
#include "Initializer/Clustering/LadderSearch.h"
#include "Initializer/Clustering/TimestepHistogram.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

namespace {

/// Hard bound on the coarsest update factor the lattice search will consider.
///
/// The dynamic program allocates O(clusters * maxIndex) state, so a badly graded mesh with a
/// huge timestep spread could otherwise ask for gigabytes. A cluster updating more than four
/// orders of magnitude less often than the finest one is not a configuration anyone wants,
/// so capping here loses nothing in practice -- cells above the cap are served by the top
/// cluster either way.
constexpr std::size_t MaxLatticeIndex = 16384;

} // namespace

SearchResult LatticeDpSearch::run(ClusteringEvaluator& evaluator,
                                  const SearchConstraints& constraints) {
  reductions_ = 0;

  if (constraints.costModel.isUpdateCount()) {
    logInfo() << "The lattice search is minimizing pure update count, which is always best "
                 "served by the finest admissible ladder. Set LtsCostLaunch and/or LtsCostFill "
                 "to make coarser ladders competitive.";
  }

  const double minimumTimestep = evaluator.timesteps().globalMinTimeStep;
  const double maximumTimestep = evaluator.timesteps().globalMaxTimeStep;

  const double minWiggleFactor = constraints.minWiggleFactor;
  const double maxWiggleFactor = 1.0;
  const double stepSizeWiggleFactor = constraints.wiggleFactorStepsize;
  const int numberOfStepsWiggleFactor =
      std::ceil((maxWiggleFactor - minWiggleFactor) / stepSizeWiggleFactor) + 1;

  double bestPredictedCost = std::numeric_limits<double>::infinity();
  double bestWiggleFactor = maxWiggleFactor;
  std::vector<std::uint64_t> bestRatios;
  bool cappedAnywhere = false;

  for (int i = 0; i < numberOfStepsWiggleFactor; ++i) {
    const double curWiggleFactor =
        std::min(minWiggleFactor + i * stepSizeWiggleFactor, maxWiggleFactor);
    const double baseTimestep = minimumTimestep * curWiggleFactor;

    const auto reachable = static_cast<std::size_t>(maximumTimestep / baseTimestep);
    const auto maxIndex = std::min(reachable, MaxLatticeIndex);
    cappedAnywhere = cappedAnywhere || maxIndex < reachable;

    const auto histogram = evaluator.timestepHistogram(curWiggleFactor, maxIndex);

    const LadderConstraints ladderConstraints{
        maxIndex, static_cast<std::size_t>(constraints.maxClusterId) + 1, constraints.maxRatio};
    const auto candidate =
        optimalLadder(histogram, constraints.costModel, baseTimestep, ladderConstraints);

    if (candidate.cost < bestPredictedCost) {
      bestPredictedCost = candidate.cost;
      bestWiggleFactor = curWiggleFactor;
      bestRatios = candidate.ratios;
    }
  }

  if (cappedAnywhere) {
    logWarning() << "The timestep spread of this mesh exceeds what the lattice search explores;"
                 << "update factors above" << MaxLatticeIndex << "were not considered.";
  }

  // Only the winner is realized: this is the expensive part, since it runs the smoothing
  // fixed point. The predicted cost ignores those demotions, so the two are logged side by
  // side -- a large gap is the signal that a shortlist would be worth re-ranking.
  reductions_ = evaluator.realize(bestRatios, bestWiggleFactor);
  const auto maxClusterId = evaluator.globalMaxClusterId();
  const auto realizedCost = evaluator.globalCost(bestRatios, bestWiggleFactor);

  logInfo() << "The best wiggle factor is" << bestWiggleFactor << "with rates" << bestRatios
            << "and" << maxClusterId + 1 << "time clusters";
  logInfo() << "Predicted cost" << bestPredictedCost << ", cost after enforcing the maximum"
            << "difference" << realizedCost;
  logInfo() << "Enforcing maximum difference for the chosen ladder took" << reductions_
            << "reductions.";

  return SearchResult{bestWiggleFactor, maxClusterId, realizedCost, bestRatios};
}

} // namespace seissol::initializer
