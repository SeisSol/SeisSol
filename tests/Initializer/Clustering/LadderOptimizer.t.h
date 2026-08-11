// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Clustering/ClusterCostModel.h"
#include "Initializer/Clustering/LadderOptimizer.h"
#include "Initializer/Clustering/TimestepHistogram.h"
#include "TestHelper.h"

#include <cstdint>
#include <mpi.h>
#include <vector>

namespace seissol::unit_test {

using seissol::initializer::ClusterCostModel;
using seissol::initializer::LadderConstraints;
using seissol::initializer::ladderCost;
using seissol::initializer::optimalLadder;
using seissol::initializer::TimestepHistogram;

using Ratios = std::vector<std::uint64_t>;

TEST_CASE("ClusterCostModel: the default is the legacy objective") {
  const ClusterCostModel model{};
  REQUIRE(model.isUpdateCount());
  // exactly weight/updateFactor, so a search configured this way optimizes what SeisSol has
  // always optimized
  REQUIRE(model.clusterTerm(1234.0, 8) == 1234.0 / 8.0);
  REQUIRE(model.clusterTerm(0.0, 4) == 0.0);

  SUBCASE("a launch cost is paid even by an empty cluster") {
    const ClusterCostModel withLaunch{1.0, 50.0, 0.0};
    REQUIRE_FALSE(withLaunch.isUpdateCount());
    REQUIRE(withLaunch.clusterTerm(0.0, 2) == 25.0);
  }

  SUBCASE("a fill threshold lifts underfilled clusters") {
    const ClusterCostModel withFill{1.0, 0.0, 400.0};
    REQUIRE(withFill.clusterTerm(100.0, 1) == 400.0);
    REQUIRE(withFill.clusterTerm(900.0, 1) == 900.0);
  }
}

TEST_CASE("TimestepHistogram: binning") {
  const auto timesteps = std::vector<double>{1.0, 1.5, 2.0, 3.9, 4.0, 100.0};
  const auto costs = std::vector<int>{1, 2, 4, 8, 16, 32};

  const auto histogram = TimestepHistogram::fromCells(timesteps, costs, 1.0, 4);
  REQUIRE(histogram.maxIndex() == 4);
  REQUIRE(histogram.totalWeight() == AbsApprox(63.0));

  REQUIRE(histogram.weightBelow(0) == AbsApprox(0.0));
  REQUIRE(histogram.weightBelow(1) == AbsApprox(0.0));
  REQUIRE(histogram.weightIn(1, 2) == AbsApprox(3.0));  // 1.0 and 1.5
  REQUIRE(histogram.weightIn(2, 3) == AbsApprox(4.0));  // 2.0
  REQUIRE(histogram.weightIn(3, 4) == AbsApprox(8.0));  // 3.9
  REQUIRE(histogram.weightIn(4, 5) == AbsApprox(48.0)); // 4.0 and the clamped 100.0

  SUBCASE("indices above the range are folded into the top bin, not wrapped") {
    const auto narrow = TimestepHistogram::fromCells(timesteps, costs, 1.0, 2);
    REQUIRE(narrow.maxIndex() == 2);
    REQUIRE(narrow.totalWeight() == AbsApprox(63.0));
    REQUIRE(narrow.weightIn(2, 3) == AbsApprox(60.0));
  }

  SUBCASE("querying past the end saturates") {
    REQUIRE(histogram.weightBelow(9999) == AbsApprox(histogram.totalWeight()));
  }

  SUBCASE("reducing over a single rank changes nothing") {
    auto reduced = histogram;
    reduced.reduce(MPI_COMM_SELF);
    REQUIRE(reduced.cumulative() == histogram.cumulative());
  }
}

TEST_CASE("optimalLadder: degenerate inputs") {
  const auto timesteps = std::vector<double>{1.0, 1.0, 1.0};
  const auto costs = std::vector<int>{5, 5, 5};
  const auto histogram = TimestepHistogram::fromCells(timesteps, costs, 1.0, 8);

  SUBCASE("everything in one bin gives a single cluster") {
    const auto result =
        optimalLadder(histogram, ClusterCostModel{}, 1.0, LadderConstraints{8, 8, 0});
    REQUIRE(result.ratios.empty());
    REQUIRE(result.cost == AbsApprox(15.0));
  }

  SUBCASE("a single admissible cluster leaves no choice") {
    const auto result =
        optimalLadder(histogram, ClusterCostModel{}, 1.0, LadderConstraints{8, 1, 0});
    REQUIRE(result.ratios.empty());
  }

  SUBCASE("maxIndex of one collapses the lattice") {
    const auto narrow = TimestepHistogram::fromCells(timesteps, costs, 1.0, 1);
    const auto result = optimalLadder(narrow, ClusterCostModel{}, 1.0, LadderConstraints{1, 8, 0});
    REQUIRE(result.ratios.empty());
  }
}

TEST_CASE("optimalLadder: splitting a bimodal distribution") {
  // half the work at the finest timestep, half at eight times that
  const auto timesteps = std::vector<double>{1.0, 1.0, 8.0, 8.0};
  const auto costs = std::vector<int>{100, 100, 100, 100};
  const auto histogram = TimestepHistogram::fromCells(timesteps, costs, 1.0, 8);

  SUBCASE("update counting takes the direct jump") {
    const auto result =
        optimalLadder(histogram, ClusterCostModel{}, 1.0, LadderConstraints{8, 8, 0});
    // 200/1 + 200/8; going via an intermediate empty rung costs the same, so the tie-break
    // towards fewer clusters has to pick the two-cluster ladder
    REQUIRE(result.ratios == Ratios{8});
    REQUIRE(result.cost == AbsApprox(225.0));
    REQUIRE(ladderCost(histogram, ClusterCostModel{}, 1.0, result.ratios) ==
            AbsApprox(result.cost));
  }

  SUBCASE("a ratio cap forces the ladder through intermediate rungs") {
    const auto result =
        optimalLadder(histogram, ClusterCostModel{}, 1.0, LadderConstraints{8, 8, 2});
    REQUIRE(result.ratios == Ratios{2, 2, 2});
    // the intermediate clusters are empty, so update counting charges nothing for them
    REQUIRE(result.cost == AbsApprox(225.0));
  }

  SUBCASE("a launch cost makes empty intermediate rungs unaffordable") {
    const ClusterCostModel model{1.0, 60.0, 0.0};
    const auto capped = optimalLadder(histogram, model, 1.0, LadderConstraints{8, 8, 2});
    const auto free = optimalLadder(histogram, model, 1.0, LadderConstraints{8, 8, 0});
    REQUIRE(free.ratios == Ratios{8});
    REQUIRE(free.cost < capped.cost);
  }

  SUBCASE("the cluster budget is respected") {
    const auto result =
        optimalLadder(histogram, ClusterCostModel{}, 1.0, LadderConstraints{8, 2, 2});
    REQUIRE(result.ratios.size() + 1 <= 2);
  }
}

TEST_CASE("optimalLadder: every ratio is a proper divisor step") {
  const auto timesteps = std::vector<double>{1.0, 2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 23.0};
  const auto costs = std::vector<int>{9, 8, 7, 6, 5, 4, 3, 2, 1};
  const auto histogram = TimestepHistogram::fromCells(timesteps, costs, 1.0, 24);

  for (const auto& model : std::vector<ClusterCostModel>{
           ClusterCostModel{}, ClusterCostModel{1.0, 3.0, 0.0}, ClusterCostModel{2.0, 1.0, 6.0}}) {
    const auto result = optimalLadder(histogram, model, 1.0, LadderConstraints{24, 8, 0});
    std::uint64_t updateFactor = 1;
    for (const auto ratio : result.ratios) {
      REQUIRE(ratio >= 2);
      updateFactor *= ratio;
    }
    REQUIRE(updateFactor <= 24);
    REQUIRE(ladderCost(histogram, model, 1.0, result.ratios) == AbsApprox(result.cost));
  }
}

} // namespace seissol::unit_test
