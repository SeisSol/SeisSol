// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERLADDER_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERLADDER_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace seissol::initializer {

/// A normalized LTS timestep ladder.
///
/// The user-facing rate vector (`ClusteredLTS` in the parameter file) is a compressed
/// description: it may be shorter than the ladder, in which case it is extended with its
/// last entry, and a ratio of 1 at an index >= 1 terminates it. This class expands that
/// description once, so that everything downstream can work with a plain list of ratios of
/// known length.
///
/// Cluster k covers cell timesteps in `[boundary(k-1), boundary(k))`, has update factor
/// `updateFactor(k) = prod_{i<k} ratios[i]` and runs at `timestep(k)`.
///
/// # Two deliberately distinct floating-point paths
///
/// `boundary()` accumulates (`b[0] = wiggle * r[0] * dtMin`, then `b[k] = b[k-1] * r[k]`)
/// while `timestep()` multiplies the exact integer update factor by the base timestep.
/// These can differ in the last ulp and are NOT interchangeable: the first reproduces the
/// binning of the legacy `getCluster()`, the second reproduces `ClusterLayout` and thereby
/// the timestep the solver actually integrates with. Unifying them would move cells across
/// cluster boundaries or perturb the simulated dt, so both are kept as they are.
class ClusterLadder {
  public:
  static constexpr std::size_t Unbounded = std::numeric_limits<std::size_t>::max();

  /// Number of clusters the rate vector itself can express, i.e. the position of the
  /// terminating ratio of 1, or `Unbounded` if there is none.
  ///
  /// Note that `rate[0]` is never a terminator -- it is an ordinary factor. `{1, 2}`
  /// therefore describes an unbounded ladder whose cluster 0 is structurally empty.
  static std::size_t intrinsicClusterCount(const std::vector<std::uint64_t>& rate);

  /// Expands `rate` into exactly `clusterCount - 1` ratios.
  static std::vector<std::uint64_t> normalize(const std::vector<std::uint64_t>& rate,
                                              std::size_t clusterCount);

  /// Ladder for binning cell timesteps. The number of clusters is the smaller of the
  /// intrinsic count and what `maximumTimestep` can reach, so that `clusterOf()` never
  /// reports a cluster that no cell can populate.
  ///
  /// `minimumTimestep` is the *pre-wiggle* global minimum.
  static ClusterLadder forBinning(const std::vector<std::uint64_t>& rate,
                                  double minimumTimestep,
                                  double wiggleFactor,
                                  double maximumTimestep);

  /// Ladder from an already complete ratio list: no repetition of the last entry, no
  /// terminator, no truncation. `ratios` *is* the ladder.
  ///
  /// This is the counterpart of `forBinning` for ladders that carry their own length, such as
  /// the output of a search. Applying `forBinning` to such a list would extend it, and that
  /// extension is not a no-op even when the list came from `forBinning` in the first place:
  /// the repeated entry is the last ratio, not the one the original vector would have used at
  /// that position.
  static ClusterLadder
      exact(const std::vector<std::uint64_t>& ratios, double minimumTimestep, double wiggleFactor);

  /// Ladder of a fixed size, for a base timestep that already includes the wiggle factor.
  /// `clusterOf()` on such a ladder is not bit-identical to the binning path (see above);
  /// it exists so that `ClusterLayout` can share the update-factor logic.
  static ClusterLadder
      ofSize(const std::vector<std::uint64_t>& rate, double baseTimestep, std::size_t clusterCount);

  [[nodiscard]] std::size_t clusterCount() const { return ratios_.size() + 1; }

  /// Normalized ratios; `ratios()[k]` is `timestep(k+1) / timestep(k)`.
  [[nodiscard]] const std::vector<std::uint64_t>& ratios() const { return ratios_; }

  /// Update factor of cluster k relative to cluster 0. Equals the legacy
  /// `ratepow(rate, 0, k)`.
  [[nodiscard]] std::uint64_t updateFactor(std::size_t cluster) const;

  /// Timestep of cluster k, wiggle factor included.
  [[nodiscard]] double timestep(std::size_t cluster) const;

  /// Upper bound of cluster k, i.e. the smallest cell timestep that lands in cluster k+1.
  /// Only defined for `k + 1 < clusterCount()`.
  [[nodiscard]] double boundary(std::size_t cluster) const { return boundaries_[cluster]; }

  /// Bins a cell timestep. Equivalent to the legacy `getCluster()`.
  [[nodiscard]] std::size_t clusterOf(double cellTimestep) const;

  [[nodiscard]] double baseTimestep() const { return baseTimestep_; }

  /// True if `ratios()[0] == 1`, in which case cluster 0 shares its timestep with cluster 1
  /// and can never be populated. Legal, and reachable via a leading 1 in the rate vector.
  [[nodiscard]] bool hasEmptyBaseCluster() const { return !ratios_.empty() && ratios_[0] == 1; }

  /// Drops the topmost clusters. Used for the cluster cap and for top merges.
  [[nodiscard]] ClusterLadder truncated(std::size_t clusterCount) const;

  private:
  ClusterLadder(std::vector<std::uint64_t> ratios,
                std::vector<double> boundaries,
                std::vector<std::uint64_t> updateFactors,
                double baseTimestep);

  std::vector<std::uint64_t> ratios_;
  std::vector<double> boundaries_;
  std::vector<std::uint64_t> updateFactors_;
  double baseTimestep_{0.0};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERLADDER_H_
