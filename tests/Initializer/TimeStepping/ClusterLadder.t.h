// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Characterization tests for the LTS ladder semantics as implemented on master.
//
// These pin the *current* behavior of getCluster()/ratepow()/ClusterLayout so that the
// upcoming ClusterLadder consolidation can be shown to be behavior-preserving. They are
// deliberately written against observable behavior, not against an intended semantics --
// where the current implementation is inconsistent, the inconsistency is pinned and
// annotated rather than "fixed" here.

#include "Initializer/TimeStepping/ClusterLadder.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/TimeStepping/LtsWeights/LtsWeights.h"
#include "TestHelper.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace seissol::unit_test {

namespace {

using Rates = std::vector<std::uint64_t>;

constexpr auto Unbounded = std::numeric_limits<std::size_t>::max();

/// Number of clusters implied by the rate-1 terminator, or Unbounded if there is none.
///
/// Mirrors the termination rule of getCluster(): at cluster c the next ratio is
/// rate[c+1], or rate.back() once c+1 runs past the end; a ratio of 1 terminates the
/// ladder *without* emitting a further cluster. Note that rate[0] is never consulted as a
/// terminator -- it is an ordinary factor -- so {1, 2} does not terminate, it merely
/// leaves cluster 0 empty.
std::size_t clusterCountOf(const Rates& rate) {
  if (rate.empty()) {
    return 1;
  }
  for (std::size_t i = 1; i < rate.size(); ++i) {
    if (rate[i] == 1) {
      return i;
    }
  }
  // the back()-extension can terminate as well; this can only bite for a single-entry {1}
  if (rate.back() == 1) {
    return rate.size();
  }
  return Unbounded;
}

/// Independent reference implementation of the binning, built from ratepow() boundaries
/// and a linear scan instead of getCluster()'s floating-point accumulation. Cross-checking
/// the two is the whole point: it ties the forward map (timestep -> cluster id) to the
/// inverse map (cluster id -> update factor) that ratepow() provides.
std::uint64_t referenceCluster(double timestep,
                               double globalMinTimestep,
                               double wiggleFactor,
                               const Rates& rate) {
  using namespace seissol::initializer::time_stepping;
  if (rate.empty()) {
    return 0;
  }
  const auto count = clusterCountOf(rate);
  const auto scaled = timestep / (wiggleFactor * globalMinTimestep);

  std::uint64_t cluster = 0;
  while (count == Unbounded || cluster + 1 < count) {
    const auto boundary = static_cast<double>(ratepow(rate, 0, cluster + 1));
    if (boundary > scaled) {
      break;
    }
    ++cluster;
  }
  return cluster;
}

std::uint64_t cluster(double timestep, const Rates& rate, double wiggle = 1.0) {
  using namespace seissol::initializer::time_stepping;
  return getCluster(timestep, 1.0, wiggle, rate);
}

/// Largest double strictly below x. Used to probe just inside a bin boundary.
double justBelow(double value) {
  return std::nextafter(value, -std::numeric_limits<double>::max());
}

} // namespace

TEST_CASE("LTS ladder: normalization table") {
  // dt_min = 1 and wiggle = 1 throughout, so every bin boundary is an exact integer in
  // double and the <= comparison in getCluster() is exact.

  SUBCASE("empty rate vector is GTS") {
    const Rates rate{};
    for (const double timestep : {0.5, 1.0, 2.0, 100.0, 1e6}) {
      REQUIRE(cluster(timestep, rate) == 0);
    }
  }

  SUBCASE("rate {1} is GTS") {
    const Rates rate{1};
    REQUIRE(clusterCountOf(rate) == 1);
    for (const double timestep : {0.5, 1.0, 2.0, 100.0, 1e6}) {
      REQUIRE(cluster(timestep, rate) == 0);
    }
  }

  SUBCASE("uniform rate {2}") {
    const Rates rate{2};
    REQUIRE(clusterCountOf(rate) == Unbounded);
    // boundaries at 2, 4, 8, 16, ...
    REQUIRE(cluster(0.5, rate) == 0);
    REQUIRE(cluster(1.0, rate) == 0);
    REQUIRE(cluster(justBelow(2.0), rate) == 0);
    REQUIRE(cluster(2.0, rate) == 1);
    REQUIRE(cluster(justBelow(4.0), rate) == 1);
    REQUIRE(cluster(4.0, rate) == 2);
    REQUIRE(cluster(1024.0, rate) == 10);
  }

  SUBCASE("bottom merge {4, 2}") {
    const Rates rate{4, 2};
    REQUIRE(clusterCountOf(rate) == Unbounded);
    // boundaries at 4, 8, 16, ... -- cluster 0 absorbs four times the usual range
    REQUIRE(cluster(justBelow(4.0), rate) == 0);
    REQUIRE(cluster(4.0, rate) == 1);
    REQUIRE(cluster(justBelow(8.0), rate) == 1);
    REQUIRE(cluster(8.0, rate) == 2);
    REQUIRE(cluster(16.0, rate) == 3);
  }

  SUBCASE("leading one leaves cluster 0 empty but does not terminate") {
    const Rates rate{1, 2};
    // rate[0] == 1 is an ordinary factor, so the ladder continues. With wiggle == 1 the
    // first bin is [0, 1) and thus unreachable for any cell (tau >= dt_min).
    REQUIRE(clusterCountOf(rate) == Unbounded);
    REQUIRE(cluster(justBelow(1.0), rate) == 0);
    REQUIRE(cluster(1.0, rate) == 1);
    REQUIRE(cluster(justBelow(2.0), rate) == 1);
    REQUIRE(cluster(2.0, rate) == 2);
  }

  SUBCASE("trailing one truncates the ladder") {
    const Rates rate{2, 2, 1};
    REQUIRE(clusterCountOf(rate) == 2);
    REQUIRE(cluster(justBelow(2.0), rate) == 0);
    REQUIRE(cluster(2.0, rate) == 1);
    // cluster 1 is the last one and absorbs everything above, even past 4
    REQUIRE(cluster(4.0, rate) == 1);
    REQUIRE(cluster(1e6, rate) == 1);
  }

  SUBCASE("mixed rates {3, 2, 5, 6, 1}") {
    const Rates rate{3, 2, 5, 6, 1};
    // Terminator at index 4 => clusters 0..3, i.e. inter-cluster ratios {3, 2, 5}.
    // rate[3] == 6 is consumed but never yields a cluster: it is a dead entry.
    REQUIRE(clusterCountOf(rate) == 4);
    REQUIRE(cluster(justBelow(3.0), rate) == 0);
    REQUIRE(cluster(3.0, rate) == 1);
    REQUIRE(cluster(justBelow(6.0), rate) == 1);
    REQUIRE(cluster(6.0, rate) == 2);
    REQUIRE(cluster(justBelow(30.0), rate) == 2);
    REQUIRE(cluster(30.0, rate) == 3);
    // top cluster is unbounded despite rate[3] == 6 suggesting a boundary at 180
    REQUIRE(cluster(180.0, rate) == 3);
    REQUIRE(cluster(1e6, rate) == 3);

    // the dead entry really is dead
    for (const std::uint64_t dead : {2, 4, 6, 7, 99}) {
      const Rates variant{3, 2, 5, dead, 1};
      for (const double timestep : {1.0, 3.0, 5.9, 6.0, 29.0, 30.0, 1e5}) {
        REQUIRE(cluster(timestep, variant) == cluster(timestep, rate));
      }
    }
  }
}

TEST_CASE("LTS ladder: getCluster and ratepow describe the same ladder") {
  // The invariant that ClusterLadder has to preserve:
  //   getCluster(tau) == k  <=>  ratepow(rate,0,k) <= tau/(w*dt) < ratepow(rate,0,k+1)
  // (with the top cluster unbounded wherever the rate-1 terminator applies).
  const std::vector<Rates> rateVectors{
      {2}, {3}, {4, 2}, {2, 3}, {1, 2}, {2, 2, 1}, {3, 2, 5, 6, 1}, {2, 4, 8}};

  for (const auto& rate : rateVectors) {
    CAPTURE(rate);
    for (const double wiggle : {1.0, 0.9, 0.51, 0.25}) {
      CAPTURE(wiggle);
      for (int i = 1; i <= 400; ++i) {
        const double timestep = 0.25 * i;
        CAPTURE(timestep);
        REQUIRE(cluster(timestep, rate, wiggle) == referenceCluster(timestep, 1.0, wiggle, rate));
      }
    }
  }
}

TEST_CASE("LTS ladder: cluster id is monotone in the wiggle factor") {
  // Underpins the clustering cache in computeClusterIdsAndEnforceMaximumDifferenceCached():
  // shrinking the wiggle factor shrinks every bin boundary, so cell cluster ids can only
  // rise. Sweeping wiggle factors upwards therefore lets min(cached, fresh) start the
  // maximum-difference fixed point from an almost-converged state.
  const std::vector<Rates> rateVectors{{2}, {4, 2}, {3, 2, 5, 6, 1}};

  for (const auto& rate : rateVectors) {
    CAPTURE(rate);
    for (int i = 1; i <= 200; ++i) {
      const double timestep = 0.5 * i;
      CAPTURE(timestep);
      std::uint64_t previous = std::numeric_limits<std::uint64_t>::max();
      // ascending wiggle => non-increasing cluster id
      for (int j = 10; j <= 100; ++j) {
        const double wiggle = 0.01 * j;
        const auto current = cluster(timestep, rate, wiggle);
        REQUIRE(current <= previous);
        previous = current;
      }
    }
  }
}

TEST_CASE("LTS ladder: ratepow and ClusterLayout::clusterRate") {
  using namespace seissol::initializer;
  using namespace seissol::initializer::time_stepping;

  SUBCASE("they agree wherever no rate equals one") {
    const std::vector<Rates> rateVectors{{2}, {3}, {4, 2}, {2, 3}, {2, 4, 8}};
    for (const auto& rate : rateVectors) {
      CAPTURE(rate);
      const ClusterLayout layout(rate, 1.0, 8);
      for (std::size_t k = 0; k < 8; ++k) {
        REQUIRE(layout.clusterRate(k) == ratepow(rate, 0, k));
        REQUIRE(layout.timestepRate(k) == AbsApprox(static_cast<double>(ratepow(rate, 0, k))));
      }
    }
  }

  SUBCASE("the terminator is now honoured on both sides") {
    // Resolved inconsistency. On master getCluster() read a ratio of 1 as a terminator
    // while ratepow()/clusterRate() read it as a literal factor, so for {2, 2, 1} the
    // ladder stopped at cluster 1 yet clusterRate(3) still reported 4 -- an update factor
    // for a cluster that can never be populated. ClusterLadder settles on getCluster()'s
    // reading: the ladder simply has no cluster 2 or 3 to ask about.
    const Rates rate{2, 2, 1};
    REQUIRE(clusterCountOf(rate) == 2);
    REQUIRE(ClusterLadder::intrinsicClusterCount(rate) == 2);

    const ClusterLayout layout(rate, 1.0, 2);
    REQUIRE(layout.globalClusterCount == 2);
    REQUIRE(layout.rates == Rates{2}); // normalized: the trailing 1 is gone
    REQUIRE(layout.clusterRate(0) == 1);
    REQUIRE(layout.clusterRate(1) == 2);
    REQUIRE(cluster(1e6, rate) == 1);
  }

  SUBCASE("both extend a short rate vector with its last entry") {
    const Rates rate{4, 2};
    const ClusterLayout layout(rate, 1.0, 6);
    REQUIRE(ratepow(rate, 0, 4) == 4 * 2 * 2 * 2);
    REQUIRE(layout.clusterRate(4) == 4 * 2 * 2 * 2);
    // partial products start at an offset, too
    REQUIRE(ratepow(rate, 1, 4) == 2 * 2 * 2);
    REQUIRE(ratepow(rate, 2, 2) == 1);
  }
}

TEST_CASE("ClusterLadder: normalization") {
  using seissol::initializer::ClusterLadder;

  SUBCASE("intrinsic cluster count") {
    REQUIRE(ClusterLadder::intrinsicClusterCount({}) == 1);
    REQUIRE(ClusterLadder::intrinsicClusterCount({1}) == 1);
    REQUIRE(ClusterLadder::intrinsicClusterCount({2, 1}) == 1);
    REQUIRE(ClusterLadder::intrinsicClusterCount({2, 2, 1}) == 2);
    REQUIRE(ClusterLadder::intrinsicClusterCount({3, 2, 5, 6, 1}) == 4);
    REQUIRE(ClusterLadder::intrinsicClusterCount({2}) == ClusterLadder::Unbounded);
    REQUIRE(ClusterLadder::intrinsicClusterCount({4, 2}) == ClusterLadder::Unbounded);
    // a leading 1 is an ordinary factor, not a terminator
    REQUIRE(ClusterLadder::intrinsicClusterCount({1, 2}) == ClusterLadder::Unbounded);
  }

  SUBCASE("short vectors are expanded with their last entry") {
    REQUIRE(ClusterLadder::normalize({4, 2}, 5) == Rates{4, 2, 2, 2});
    REQUIRE(ClusterLadder::normalize({2}, 4) == Rates{2, 2, 2});
    REQUIRE(ClusterLadder::normalize({3, 2, 5, 6, 1}, 4) == Rates{3, 2, 5});
    REQUIRE(ClusterLadder::normalize({2}, 1).empty());
  }

  SUBCASE("a leading one survives normalization and marks an empty base cluster") {
    const auto ladder = ClusterLadder::forBinning({1, 2}, 1.0, 1.0, 100.0);
    REQUIRE(ladder.hasEmptyBaseCluster());
    REQUIRE(ladder.ratios()[0] == 1);
    REQUIRE(ladder.timestep(0) == AbsApprox(ladder.timestep(1)));
    // no cell can reach cluster 0, since every cell timestep is at least dt_min
    REQUIRE(ladder.clusterOf(1.0) == 1);
  }
}

TEST_CASE("ClusterLadder: agrees with the legacy free functions") {
  using namespace seissol::initializer::time_stepping;
  using seissol::initializer::ClusterLadder;

  const std::vector<Rates> rateVectors{
      {2}, {3}, {4, 2}, {2, 3}, {1, 2}, {2, 2, 1}, {3, 2, 5, 6, 1}, {2, 4, 8}};
  constexpr double MaximumTimestep = 750.0;

  for (const auto& rate : rateVectors) {
    CAPTURE(rate);
    for (const double wiggle : {1.0, 0.9, 0.51, 0.25}) {
      CAPTURE(wiggle);
      const auto ladder = ClusterLadder::forBinning(rate, 1.0, wiggle, MaximumTimestep);

      // the ladder stops exactly where the largest cell timestep does
      REQUIRE(ladder.clusterCount() == getCluster(MaximumTimestep, 1.0, wiggle, rate) + 1);

      for (std::size_t k = 0; k < ladder.clusterCount(); ++k) {
        REQUIRE(ladder.updateFactor(k) == ratepow(rate, 0, k));
      }
      // beyond index 0 the normalized ratios never contain a 1
      for (std::size_t k = 1; k < ladder.ratios().size(); ++k) {
        REQUIRE(ladder.ratios()[k] != 1);
      }

      for (int i = 1; i <= 3000; ++i) {
        const double timestep = 0.25 * i;
        CAPTURE(timestep);
        REQUIRE(ladder.clusterOf(timestep) == getCluster(timestep, 1.0, wiggle, rate));
      }
    }
  }
}

TEST_CASE("ClusterLadder: truncation") {
  using seissol::initializer::ClusterLadder;
  const auto ladder = ClusterLadder::forBinning({4, 2}, 1.0, 1.0, 1000.0);
  REQUIRE(ladder.clusterCount() > 3);

  const auto truncated = ladder.truncated(3);
  REQUIRE(truncated.clusterCount() == 3);
  REQUIRE(truncated.ratios() == Rates{4, 2});
  for (std::size_t k = 0; k < 3; ++k) {
    REQUIRE(truncated.updateFactor(k) == ladder.updateFactor(k));
    REQUIRE(truncated.timestep(k) == AbsApprox(ladder.timestep(k)));
  }
  // everything above the new top cluster collapses into it
  REQUIRE(truncated.clusterOf(1e6) == 2);
  // truncating to at least the current size is a no-op
  REQUIRE(ladder.truncated(ladder.clusterCount() + 5).clusterCount() == ladder.clusterCount());
}

} // namespace seissol::unit_test
