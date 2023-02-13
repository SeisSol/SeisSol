#include "Geometry/PUMLReader.h"
#include "Initializer/time_stepping/LtsWeights/WeightsModels.h"
#include <memory>
#include <numeric>

namespace seissol::unit_test {

TEST_CASE("LTS Weights") {
// PUMLReader is only available with MPI
#ifdef USE_MPI
  std::cout.setstate(std::ios_base::failbit);
  using namespace seissol::initializers::time_stepping;
  LtsWeightsConfig config{"Testing/material.yaml", 2, 1, 1, 1};

  auto ltsParameters = std::make_unique<LtsParameters>(2, 1.0, 0.01, false, 100, false, 1.0);
  auto ltsWeights = std::make_unique<ExponentialWeights>(config, ltsParameters.get());
  seissol::PUMLReader pumlReader("Testing/mesh.h5", 5000.0, "", ltsWeights.get());
  std::cout.clear();

  std::array<unsigned, 24> expectedWeights = {2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2,
                                              2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1};
  for (int i = 0; i < 24; i++) {
    REQUIRE(ltsWeights->vertexWeights()[i] == expectedWeights[i]);
  }
#endif
}

TEST_CASE("Cost function for LTS") {
  const auto eps = 10e-12;
  using namespace initializers::time_stepping;

  SUBCASE("No clusters") {
    std::vector<int> clusterIds = {};
    std::vector<int> cellCosts = {};
    const auto is = computeLocalCostOfClustering(clusterIds, cellCosts, 2, 1.0, 1.0);
    const auto should = 0.0;
    REQUIRE(AbsApprox(is).epsilon(eps) == should);
  }

  SUBCASE("One cluster") {
    std::vector<int> clusterIds = {0, 0, 0, 0, 0};
    std::vector<int> cellCosts = {1, 2, 3, 4, 5};
    for (int i = 1; i <= 10; ++i) {
      const auto dt = 1.0 / i;
      for (int j = 1; j <= 10; ++j) {
        const auto wiggleFactor = 1.0 / j;

        const auto is = computeLocalCostOfClustering(clusterIds, cellCosts, 2, wiggleFactor, dt);

        const auto totalCost = std::accumulate(cellCosts.begin(), cellCosts.end(), 0);
        const auto effectiveDt = dt * wiggleFactor;

        const auto should = totalCost * (1.0 / effectiveDt);
        REQUIRE(AbsApprox(is).epsilon(eps) == should);
      }
    }
  }

  SUBCASE("Two clusters") {
    std::vector<int> clusterIds = {1, 0, 1, 1};
    std::vector<int> cellCosts = {2, 1, 3, 1};
    const auto cellCostsCluster0 = 1;
    const auto cellCostsCluster1 = 2 + 3 + 1;
    for (unsigned int rate = 1; rate < 4; ++rate) {
      for (int i = 1; i <= 10; ++i) {
        const auto dt = 1.0 / i;
        for (int j = 1; j <= 10; ++j) {
          const auto wiggleFactor = 1.0 / j;

          const auto is =
              computeLocalCostOfClustering(clusterIds, cellCosts, rate, wiggleFactor, dt);

          const auto effectiveDtCluster0 = dt * wiggleFactor;
          const auto effectiveDtCluster1 = rate * effectiveDtCluster0;

          const auto costCluster0 = cellCostsCluster0 * (1.0 / effectiveDtCluster0);
          const auto costCluster1 = cellCostsCluster1 * (1.0 / effectiveDtCluster1);
          const auto should = costCluster0 + costCluster1;
          REQUIRE(AbsApprox(is).epsilon(eps) == should);
        }
      }
    }
  }

  SUBCASE("Three clusters") {
    std::vector<int> clusterIds = {2, 0, 1, 1, 1};
    std::vector<int> cellCosts = {2, 1, 3, 1, 2};
    const auto cellCostsCluster0 = 1;
    const auto cellCostsCluster1 = 2 + 3 + 1;
    const auto cellCostsCluster2 = 2;
    for (unsigned int rate = 1; rate < 4; ++rate) {
      for (int i = 1; i <= 10; ++i) {
        const auto dt = 1.0 / i;
        for (int j = 1; j <= 10; ++j) {
          const auto wiggleFactor = 1.0 / j;

          const auto is =
              computeLocalCostOfClustering(clusterIds, cellCosts, rate, wiggleFactor, dt);

          const auto effectiveDtCluster0 = dt * wiggleFactor;
          const auto effectiveDtCluster1 = rate * effectiveDtCluster0;
          const auto effectiveDtCluster2 = rate * effectiveDtCluster1;

          const auto costCluster0 = cellCostsCluster0 * (1.0 / effectiveDtCluster0);
          const auto costCluster1 = cellCostsCluster1 * (1.0 / effectiveDtCluster1);
          const auto costCluster2 = cellCostsCluster2 * (1.0 / effectiveDtCluster2);
          const auto should = costCluster0 + costCluster1 + costCluster2;
          REQUIRE(AbsApprox(is).epsilon(eps) == should);
        }
      }
    }
  }
}

TEST_CASE("Enforce max cluster id") {
  using namespace seissol::initializers::time_stepping;
  const auto clusterIds = std::vector<int>{0, 1, 2, 3, 4, 5, 6, 6, 5, 4, 3, 2, 1, 0};
  SUBCASE("No change") {
    const auto should = clusterIds;
    const auto is = enforceMaxClusterId(clusterIds, 6);
    REQUIRE(is == should);
  }
  SUBCASE("Only one cluster") {
    const auto should = std::vector<int>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const auto is = enforceMaxClusterId(clusterIds, 0);
    REQUIRE(is == should);
  }
  SUBCASE("Three clusters") {
    const auto should = std::vector<int>{0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0};
    const auto is = enforceMaxClusterId(clusterIds, 2);
    REQUIRE(is == should);
  }
}

TEST_CASE("Auto merging of clusters") {
  using namespace seissol::initializers::time_stepping;
  const auto clusterIds = std::vector<int>{0, 0, 0, 0, 1, 1, 2};
  const auto cellCosts = std::vector<int>{1, 1, 1, 1, 3, 3, 9};
  const auto minDt = 0.5;
  const auto costBefore = computeLocalCostOfClustering(clusterIds, cellCosts, 2, 1.0, minDt);

  SUBCASE("Reduces to GTS") {
    const auto should = 0;
    const auto is = computeMaxClusterIdAfterAutoMerge(
        clusterIds, cellCosts, 2, std::numeric_limits<double>::max(), 1.0, minDt);
    REQUIRE(is == should);
  }

  SUBCASE("Does nothing for GTS") {
    const auto should = 0;
    for (int i = 1; i <= 5; ++i) {
      const auto is = computeMaxClusterIdAfterAutoMerge(
          enforceMaxClusterId(clusterIds, 0), cellCosts, 1, i, 0, 0);
      REQUIRE(is == should);
    }
  }

  SUBCASE("No performance loss allowed") {
    const auto should = 2;
    SUBCASE("Rate 2") {
      const auto is = computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, 2, costBefore, 1, minDt);
      REQUIRE(is == should);
    }
    SUBCASE("Rate 3") {
      const auto is = computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, 3, costBefore, 1, minDt);
      REQUIRE(is == should);
    }
  }

  SUBCASE("Some performance loss allowed") {
    SUBCASE("Merge one cluster") {
      const auto should = 1;
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, 2, 1.25 * costBefore, 1, minDt);
      REQUIRE(is == should);
    }
    SUBCASE("Merge two clusters") {
      const auto should = 0;
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, 2, 2.06 * costBefore, 1, minDt);
      REQUIRE(is == should);
    }
  }
}

} // namespace seissol::unit_test