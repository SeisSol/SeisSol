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

  auto ltsWeights = std::make_unique<ExponentialWeights>(config);
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
    const auto is = computeCostOfClustering(clusterIds, cellCosts, 2, 1.0, 1.0);
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

        const auto is = computeCostOfClustering(clusterIds, cellCosts, 2, wiggleFactor, dt);

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

          const auto is = computeCostOfClustering(clusterIds, cellCosts, rate, wiggleFactor, dt);

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

          const auto is = computeCostOfClustering(clusterIds, cellCosts, rate, wiggleFactor, dt);

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

} // namespace seissol::unit_test