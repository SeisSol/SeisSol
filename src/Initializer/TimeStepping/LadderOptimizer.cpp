// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/TimeStepping/LadderOptimizer.h"

#include "Initializer/TimeStepping/ClusterCostModel.h"
#include "Initializer/TimeStepping/TimestepHistogram.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace seissol::initializer {

namespace {

/// Longest divisibility chain that fits below `maxIndex`. Ratios are at least two, so the
/// chain doubles at every step.
std::size_t longestChain(std::size_t maxIndex) {
  std::size_t rungs = 1;
  for (std::size_t factor = 1; factor <= maxIndex / 2; factor *= 2) {
    ++rungs;
  }
  return rungs;
}

/// Largest multiple of `rung` that a transition may reach.
std::size_t transitionLimit(std::size_t rung, const LadderConstraints& constraints) {
  if (constraints.maxRatio == 0) {
    return constraints.maxIndex;
  }
  // guard the multiplication: maxRatio comes from user input
  if (rung > constraints.maxIndex / constraints.maxRatio) {
    return constraints.maxIndex;
  }
  return std::min<std::size_t>(constraints.maxIndex, rung * constraints.maxRatio);
}

} // namespace

double ladderCost(const TimestepHistogram& histogram,
                  const ClusterCostModel& costModel,
                  double baseTimestep,
                  const std::vector<std::uint64_t>& ratios) {
  assert(baseTimestep > 0.0);

  double cost = 0.0;
  std::uint64_t updateFactor = 1;
  for (const auto ratio : ratios) {
    const auto next = updateFactor * ratio;
    cost += costModel.clusterTerm(histogram.weightIn(updateFactor, next), updateFactor);
    updateFactor = next;
  }
  cost += costModel.clusterTerm(histogram.totalWeight() - histogram.weightBelow(updateFactor),
                                updateFactor);
  return cost / baseTimestep;
}

LadderCandidate optimalLadder(const TimestepHistogram& histogram,
                              const ClusterCostModel& costModel,
                              double baseTimestep,
                              const LadderConstraints& constraints) {
  assert(baseTimestep > 0.0);

  const auto maxIndex = std::max<std::size_t>(constraints.maxIndex, 1);
  const auto maxRungs = std::max<std::size_t>(
      1, std::min(longestChain(maxIndex), std::max<std::size_t>(constraints.maxClusterCount, 1)));

  const auto total = histogram.totalWeight();
  const auto stride = maxIndex + 1;

  // best[rungs][rung]: cheapest way to cover [rung, inf) with at most `rungs` clusters, the
  // finest of which has update factor `rung`
  std::vector<double> best((maxRungs + 1) * stride, std::numeric_limits<double>::infinity());
  std::vector<std::size_t> successor((maxRungs + 1) * stride, 0);
  // number of clusters the optimal solution at a state uses, so that ties can be broken
  // towards the shorter ladder. Preferring the terminal over a transition is not enough:
  // among equally priced transitions the finest step would otherwise win and pad the ladder
  // with clusters that carry nothing.
  std::vector<std::size_t> clusterCount((maxRungs + 1) * stride, 1);

  auto at = [stride](std::size_t rungs, std::size_t rung) { return rungs * stride + rung; };

  auto effectiveConstraints = constraints;
  effectiveConstraints.maxIndex = maxIndex;

  for (std::size_t rungs = 1; rungs <= maxRungs; ++rungs) {
    for (std::size_t rung = maxIndex; rung >= 1; --rung) {
      // terminal: this rung carries everything above it. Evaluated first so that ties keep
      // the shorter chain.
      auto bestCost = costModel.clusterTerm(total - histogram.weightBelow(rung), rung);
      std::size_t bestSuccessor = 0;
      std::size_t bestCount = 1;

      if (rungs > 1) {
        const auto limit = transitionLimit(rung, effectiveConstraints);
        for (std::size_t next = 2 * rung; next <= limit; next += rung) {
          const auto candidate = costModel.clusterTerm(histogram.weightIn(rung, next), rung) +
                                 best[at(rungs - 1, next)];
          const auto candidateCount = 1 + clusterCount[at(rungs - 1, next)];
          if (candidate < bestCost || (candidate == bestCost && candidateCount < bestCount)) {
            bestCost = candidate;
            bestSuccessor = next;
            bestCount = candidateCount;
          }
        }
      }

      best[at(rungs, rung)] = bestCost;
      successor[at(rungs, rung)] = bestSuccessor;
      clusterCount[at(rungs, rung)] = bestCount;
    }
  }

  std::vector<std::uint64_t> ratios;
  std::size_t rung = 1;
  for (auto rungs = maxRungs; rungs >= 1; --rungs) {
    const auto next = successor[at(rungs, rung)];
    if (next == 0) {
      break;
    }
    ratios.push_back(static_cast<std::uint64_t>(next / rung));
    rung = next;
  }

  return LadderCandidate{std::move(ratios), best[at(maxRungs, 1)] / baseTimestep};
}

} // namespace seissol::initializer
