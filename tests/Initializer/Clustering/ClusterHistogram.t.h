// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Clustering/ClusterHistogram.h"
#include "Initializer/Clustering/ClusterLadder.h"
#include "Initializer/Clustering/ClusteringCost.h"
#include "TestHelper.h"

#include <cstdint>
#include <mpi.h>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("ClusterHistogram: binning") {
  using seissol::initializer::ClusterHistogram;

  const auto clusterIds = std::vector<int>{0, 0, 0, 0, 1, 1, 2};
  const auto cellCosts = std::vector<int>{1, 1, 1, 1, 3, 3, 9};

  const auto histogram = ClusterHistogram::fromClustering(clusterIds, cellCosts, 3);
  REQUIRE(histogram.clusterCount() == 3);
  REQUIRE(histogram.weight(0) == AbsApprox(4.0));
  REQUIRE(histogram.weight(1) == AbsApprox(6.0));
  REQUIRE(histogram.weight(2) == AbsApprox(9.0));
  REQUIRE(histogram.totalWeight() == AbsApprox(19.0));

  SUBCASE("trailing clusters may be empty") {
    const auto padded = ClusterHistogram::fromClustering(clusterIds, cellCosts, 5);
    REQUIRE(padded.clusterCount() == 5);
    REQUIRE(padded.weight(3) == AbsApprox(0.0));
    REQUIRE(padded.weight(4) == AbsApprox(0.0));
    REQUIRE(padded.totalWeight() == AbsApprox(histogram.totalWeight()));
  }

  SUBCASE("reducing over a single rank changes nothing") {
    auto reduced = histogram;
    reduced.reduce(MPI_COMM_SELF);
    REQUIRE(reduced.weights() == histogram.weights());
  }
}

TEST_CASE("ClusterHistogram: cost agrees with the per-cell cost") {
  using namespace seissol::initializer;
  using seissol::initializer::ClusterHistogram;
  using seissol::initializer::ClusterLadder;

  const auto clusterIds = std::vector<int>{0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 4};
  const auto cellCosts = std::vector<int>{7, 3, 5, 1, 3, 3, 9, 2, 4, 6, 8, 5};
  constexpr std::size_t ClusterCount = 5;
  constexpr double MinimalTimestep = 4.2e-4;

  const auto histogram = ClusterHistogram::fromClustering(clusterIds, cellCosts, ClusterCount);

  SUBCASE("bit-identical for power-of-two rates") {
    // every term is a dyadic rational, so summing per cluster and summing per cell round
    // exactly the same way
    for (const auto& rate : std::vector<std::vector<std::uint64_t>>{{2}, {4}, {4, 2}}) {
      CAPTURE(rate);
      for (const double wiggle : {1.0, 0.75, 0.5}) {
        const auto ladder =
            ClusterLadder::ofSize(rate, MinimalTimestep * wiggle, histogram.clusterCount());
        for (std::size_t cap = 0; cap < ClusterCount; ++cap) {
          const auto expected =
              computeLocalCostOfClustering(enforceMaxClusterId(clusterIds, static_cast<int>(cap)),
                                           cellCosts,
                                           rate,
                                           wiggle,
                                           MinimalTimestep);
          REQUIRE(histogram.cost(ladder, cap) == expected);
        }
      }
    }
  }

  SUBCASE("within a few ulp for other rates") {
    // documented deviation: the per-cell sum rounds every term separately
    for (const auto& rate : std::vector<std::vector<std::uint64_t>>{{3}, {5}, {3, 2}}) {
      CAPTURE(rate);
      const auto ladder = ClusterLadder::ofSize(rate, MinimalTimestep, histogram.clusterCount());
      for (std::size_t cap = 0; cap < ClusterCount; ++cap) {
        const auto expected =
            computeLocalCostOfClustering(enforceMaxClusterId(clusterIds, static_cast<int>(cap)),
                                         cellCosts,
                                         rate,
                                         1.0,
                                         MinimalTimestep);
        REQUIRE(histogram.cost(ladder, cap) == AbsApprox(expected).epsilon(0.0).delta(1e-12));
      }
    }
  }

  SUBCASE("an uncapped cost is the cost at the topmost cluster") {
    const auto rate = std::vector<std::uint64_t>{2};
    const auto ladder = ClusterLadder::ofSize(rate, MinimalTimestep, histogram.clusterCount());
    REQUIRE(histogram.cost(ladder) == histogram.cost(ladder, ClusterCount - 1));
  }
}

} // namespace seissol::unit_test
