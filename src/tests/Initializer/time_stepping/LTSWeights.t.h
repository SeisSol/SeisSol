// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Parameters/LtsParameters.h"
#include "tests/TestHelper.h"
#include <memory>
#include <numeric>

#include "Geometry/PUMLReader.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/TimeStepping/LtsWeights/WeightsModels.h"
#include "Initializer/Typedefs.h"
#include "SeisSol.h"

namespace seissol::unit_test {

TEST_CASE("LTS Weights") {
// PUMLReader is only available with MPI
#ifdef USE_MPI
  std::cout.setstate(std::ios_base::failbit);
  using namespace seissol::initializer::time_stepping;
  const LtsWeightsConfig config{seissol::initializer::parameters::BoundaryFormat::I32, 2, 1, 1, 1};

  const seissol::initializer::parameters::LtsParameters ltsParameters(
      2,
      1.0,
      0.01,
      false,
      100,
      false,
      1.0,
      seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor,
      seissol::initializer::parameters::LtsWeightsTypes::ExponentialWeights);
  seissol::initializer::parameters::SeisSolParameters seissolParameters;
  seissolParameters.timeStepping.lts = ltsParameters;
  seissolParameters.timeStepping.cfl = 1;
  seissolParameters.timeStepping.maxTimestepWidth = 5000.0;
  seissolParameters.model.materialFileName = "Testing/material.yaml";
  seissol::SeisSol seissolInstance(seissolParameters);

  auto ltsWeights = std::make_unique<ExponentialWeights>(config, seissolInstance);
  const auto pumlReader =
      seissol::geometry::PUMLReader("Testing/mesh.h5",
                                    "Default",
                                    seissol::initializer::parameters::BoundaryFormat::I32,
                                    ltsWeights.get());
  std::cout.clear();

  const auto givenWeights =
      std::vector<unsigned>(ltsWeights->vertexWeights(), ltsWeights->vertexWeights() + 24);

  const auto expectedWeights =
      std::vector<unsigned>{2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1};

  REQUIRE(givenWeights == expectedWeights);
#endif
}

TEST_CASE("Cost function for LTS") {
  const auto eps = 10e-12;
  using namespace initializer::time_stepping;

  SUBCASE("No clusters") {
    const std::vector<int> clusterIds = {};
    const std::vector<int> cellCosts = {};
    const auto is = computeLocalCostOfClustering(clusterIds, cellCosts, 2, 1.0, 1.0);
    const auto should = 0.0;
    REQUIRE(AbsApprox(is).epsilon(eps) == should);
  }

  SUBCASE("One cluster") {
    const std::vector<int> clusterIds = {0, 0, 0, 0, 0};
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
    const std::vector<int> clusterIds = {1, 0, 1, 1};
    const std::vector<int> cellCosts = {2, 1, 3, 1};
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
    const std::vector<int> clusterIds = {2, 0, 1, 1, 1};
    const std::vector<int> cellCosts = {2, 1, 3, 1, 2};
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
  using namespace seissol::initializer::time_stepping;
  const auto clusterIds = std::vector<int>{0, 1, 2, 3, 4, 5, 6, 6, 5, 4, 3, 2, 1, 0};
  SUBCASE("No change") {
    const auto& should = clusterIds;
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
  using namespace seissol::initializer::time_stepping;
  const auto clusterIds = std::vector<int>{0, 0, 0, 0, 1, 1, 2};
  const auto cellCosts = std::vector<int>{1, 1, 1, 1, 3, 3, 9};
  const auto minDt = 0.5;
  const auto costBeforeRate2 = computeLocalCostOfClustering(clusterIds, cellCosts, 2, 1.0, minDt);
  const auto costBeforeRate3 = computeLocalCostOfClustering(clusterIds, cellCosts, 3, 1.0, minDt);

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
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, 2, costBeforeRate2, 1, minDt);
      REQUIRE(is == should);
    }
    SUBCASE("Rate 3") {
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, 3, costBeforeRate3, 1, minDt);
      REQUIRE(is == should);
    }
  }

  SUBCASE("Some performance loss allowed") {
    SUBCASE("Merge one cluster") {
      const auto should = 1;
      const auto is = computeMaxClusterIdAfterAutoMerge(
          clusterIds, cellCosts, 2, 1.25 * costBeforeRate2, 1, minDt);
      REQUIRE(is == should);
    }
    SUBCASE("Merge two clusters") {
      const auto should = 0;
      const auto is = computeMaxClusterIdAfterAutoMerge(
          clusterIds, cellCosts, 2, 2.06 * costBeforeRate2, 1, minDt);
      REQUIRE(is == should);
    }
  }
}

} // namespace seissol::unit_test
