// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/TimeStepping/ClusterLadder.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace seissol::initializer {

namespace {

/// The ratio that takes cluster k to cluster k+1, with the last entry repeated.
std::uint64_t ratioAt(const std::vector<std::uint64_t>& rate, std::size_t cluster) {
  return cluster < rate.size() ? rate[cluster] : rate.back();
}

/// The one accumulation of bin boundaries. The association order reproduces the legacy
/// `getCluster()` exactly: `wiggle * r[0]` first, then `* dtMin`, then repeated `*= r[k]`.
double advance(double boundary, std::uint64_t ratio) {
  return boundary * static_cast<double>(ratio);
}

double firstBoundary(std::uint64_t ratio, double minimumTimestep, double wiggleFactor) {
  return wiggleFactor * static_cast<double>(ratio) * minimumTimestep;
}

} // namespace

std::size_t ClusterLadder::intrinsicClusterCount(const std::vector<std::uint64_t>& rate) {
  if (rate.empty()) {
    return 1;
  }
  // rate[0] is an ordinary factor and never terminates
  for (std::size_t i = 1; i < rate.size(); ++i) {
    if (rate[i] == 1) {
      return i;
    }
  }
  // the repeated last entry can terminate too; only reachable for a single-entry {1}
  if (rate.back() == 1) {
    return rate.size();
  }
  return Unbounded;
}

std::vector<std::uint64_t> ClusterLadder::normalize(const std::vector<std::uint64_t>& rate,
                                                    std::size_t clusterCount) {
  std::vector<std::uint64_t> ratios;
  if (rate.empty() || clusterCount <= 1) {
    return ratios;
  }
  assert(clusterCount != Unbounded);
  assert(clusterCount <= intrinsicClusterCount(rate));
  ratios.reserve(clusterCount - 1);
  for (std::size_t k = 0; k + 1 < clusterCount; ++k) {
    ratios.push_back(ratioAt(rate, k));
  }
  return ratios;
}

ClusterLadder::ClusterLadder(std::vector<std::uint64_t> ratios,
                             std::vector<double> boundaries,
                             std::vector<std::uint64_t> updateFactors,
                             double baseTimestep)
    : ratios_(std::move(ratios)), boundaries_(std::move(boundaries)),
      updateFactors_(std::move(updateFactors)), baseTimestep_(baseTimestep) {
  assert(boundaries_.size() == ratios_.size());
  assert(updateFactors_.size() == ratios_.size() + 1);
}

namespace {

std::vector<std::uint64_t> buildUpdateFactors(const std::vector<std::uint64_t>& ratios) {
  std::vector<std::uint64_t> factors;
  factors.reserve(ratios.size() + 1);
  std::uint64_t factor = 1;
  factors.push_back(factor);
  for (const auto ratio : ratios) {
    factor *= ratio;
    factors.push_back(factor);
  }
  return factors;
}

} // namespace

ClusterLadder ClusterLadder::forBinning(const std::vector<std::uint64_t>& rate,
                                        double minimumTimestep,
                                        double wiggleFactor,
                                        double maximumTimestep) {
  const auto intrinsic = intrinsicClusterCount(rate);

  std::vector<std::uint64_t> ratios;
  std::vector<double> boundaries;

  if (!rate.empty()) {
    // Grow the ladder while the topmost boundary is still reachable by some cell. Stopping
    // here is what makes the ladder finite for an unbounded rate vector, and it is a no-op
    // for the binning itself: no cell timestep exceeds `maximumTimestep`.
    for (std::size_t k = 0; intrinsic == Unbounded || k + 1 < intrinsic; ++k) {
      const auto ratio = ratioAt(rate, k);
      const auto boundary = boundaries.empty() ? firstBoundary(ratio, minimumTimestep, wiggleFactor)
                                               : advance(boundaries.back(), ratio);
      // A boundary above the largest cell timestep separates nothing, so the cluster it
      // would open stays empty. Dropping it is what makes `clusterCount()` agree with
      // `getCluster(maximumTimestep) + 1`.
      if (boundary > maximumTimestep) {
        break;
      }
      ratios.push_back(ratio);
      boundaries.push_back(boundary);
    }
  }

  auto updateFactors = buildUpdateFactors(ratios);
  return ClusterLadder(std::move(ratios),
                       std::move(boundaries),
                       std::move(updateFactors),
                       minimumTimestep * wiggleFactor);
}

ClusterLadder ClusterLadder::exact(const std::vector<std::uint64_t>& ratios,
                                   double minimumTimestep,
                                   double wiggleFactor) {
  std::vector<double> boundaries;
  boundaries.reserve(ratios.size());
  for (std::size_t k = 0; k < ratios.size(); ++k) {
    // same association order as forBinning, so that a ladder rebuilt from its own ratios bins
    // every cell exactly as the original did
    boundaries.push_back(k == 0 ? firstBoundary(ratios[0], minimumTimestep, wiggleFactor)
                                : advance(boundaries.back(), ratios[k]));
  }

  auto updateFactors = buildUpdateFactors(ratios);
  return ClusterLadder(
      ratios, std::move(boundaries), std::move(updateFactors), minimumTimestep * wiggleFactor);
}

ClusterLadder ClusterLadder::ofSize(const std::vector<std::uint64_t>& rate,
                                    double baseTimestep,
                                    std::size_t clusterCount) {
  auto ratios = normalize(rate, clusterCount);

  std::vector<double> boundaries;
  boundaries.reserve(ratios.size());
  for (std::size_t k = 0; k < ratios.size(); ++k) {
    boundaries.push_back(k == 0 ? firstBoundary(ratios[0], baseTimestep, 1.0)
                                : advance(boundaries.back(), ratios[k]));
  }

  auto updateFactors = buildUpdateFactors(ratios);
  return ClusterLadder(
      std::move(ratios), std::move(boundaries), std::move(updateFactors), baseTimestep);
}

std::uint64_t ClusterLadder::updateFactor(std::size_t cluster) const {
  assert(cluster < updateFactors_.size());
  return updateFactors_[cluster];
}

double ClusterLadder::timestep(std::size_t cluster) const {
  return static_cast<double>(updateFactor(cluster)) * baseTimestep_;
}

std::size_t ClusterLadder::clusterOf(double cellTimestep) const {
  // linear scan, in the same comparison order as the legacy loop; ladders are short
  std::size_t cluster = 0;
  while (cluster < boundaries_.size() && boundaries_[cluster] <= cellTimestep) {
    ++cluster;
  }
  return cluster;
}

ClusterLadder ClusterLadder::truncated(std::size_t clusterCount) const {
  if (clusterCount >= this->clusterCount()) {
    return *this;
  }
  assert(clusterCount >= 1);
  const auto ratioCount = clusterCount - 1;
  return ClusterLadder(
      std::vector<std::uint64_t>(ratios_.begin(), ratios_.begin() + ratioCount),
      std::vector<double>(boundaries_.begin(), boundaries_.begin() + ratioCount),
      std::vector<std::uint64_t>(updateFactors_.begin(), updateFactors_.begin() + clusterCount),
      baseTimestep_);
}

} // namespace seissol::initializer
