// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "Geometry/PUMLReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/FaceMap.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/TimeStepping/ClusterLadder.h"
#include "Initializer/TimeStepping/LtsWeights/WeightsModels.h"
#include "Initializer/Typedefs.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "TestHelper.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("LTS Weights") {
  std::cout.setstate(std::ios_base::failbit);
  using namespace seissol::initializer::time_stepping;
  const LtsWeightsConfig config{
      seissol::initializer::parameters::BoundaryFormat::I32, {2}, 1, 1, 1};

  const seissol::initializer::parameters::LtsParameters ltsParameters(
      {2},
      1.0,
      0.01,
      false,
      100,
      false,
      1.0,
      seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor,
      seissol::initializer::parameters::LtsWeightsTypes::ExponentialWeights);
  seissol::initializer::parameters::SeisSolParameters seissolParameters{};
  seissolParameters.timeStepping.lts = ltsParameters;
  seissolParameters.timeStepping.cfl = 1;
  seissolParameters.timeStepping.maxTimestepWidth = 5000.0;
  seissolParameters.model.materialFileName = tpath("Testing/material.yaml");
  seissolParameters.model.useCellHomogenizedMaterial = false;
  seissolParameters.model.plasticity = false;
  const utils::Env env("SEISSOL_");
  seissol::SeisSol seissolInstance(seissolParameters, env);

  auto ltsWeights = std::make_unique<ExponentialWeights>(config, seissolInstance);
  const auto faceMap = defaultFaceMap();
  const auto pumlReader =
      seissol::geometry::PUMLReader(tpath("Testing/mesh.h5"),
                                    "Default",
                                    faceMap,
                                    seissol::initializer::parameters::BoundaryFormat::I32,
                                    seissol::initializer::parameters::TopologyFormat::Geometric,
                                    ltsWeights.get());
  std::cout.clear();

  const auto givenWeights =
      std::vector<unsigned>(ltsWeights->vertexWeights(), ltsWeights->vertexWeights() + 24);

  const auto expectedWeights =
      std::vector<unsigned>{2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1};

  REQUIRE(givenWeights == expectedWeights);
}

TEST_CASE("Cost function for LTS") {
  const auto eps = 10e-12;
  using namespace initializer::time_stepping;

  SUBCASE("No clusters") {
    const std::vector<int> clusterIds = {};
    const std::vector<int> cellCosts = {};
    const auto is = computeLocalCostOfClustering(clusterIds, cellCosts, {2}, 1.0, 1.0);
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

        const auto is = computeLocalCostOfClustering(clusterIds, cellCosts, {2}, wiggleFactor, dt);

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
              computeLocalCostOfClustering(clusterIds, cellCosts, {rate}, wiggleFactor, dt);

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
              computeLocalCostOfClustering(clusterIds, cellCosts, {rate}, wiggleFactor, dt);

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

TEST_CASE("Batched costs of capped clusterings") {
  using namespace seissol::initializer::time_stepping;
  const auto clusterIds = std::vector<int>{0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 4};
  const auto cellCosts = std::vector<int>{7, 3, 5, 1, 3, 3, 9, 2, 4, 6, 8, 5};
  const auto maxClusterId = 4;

  // The batched form has to reproduce the one-cap-at-a-time form exactly, not just closely:
  // the resulting costs are compared against an admissibility threshold, so a last-ulp
  // difference could flip a merge decision.
  for (const auto& rate :
       std::vector<std::vector<std::uint64_t>>{{2}, {3}, {4, 2}, {3, 2, 5, 6, 1}}) {
    CAPTURE(rate);
    for (const double wiggle : {1.0, 0.87, 0.5}) {
      CAPTURE(wiggle);
      for (const double minimalTimestep : {1.0, 4.2e-4}) {
        const auto batched = computeLocalCostsOfCappedClusterings(
            clusterIds, cellCosts, rate, wiggle, minimalTimestep, maxClusterId);
        REQUIRE(batched.size() == static_cast<std::size_t>(maxClusterId) + 1);
        for (int cap = 0; cap <= maxClusterId; ++cap) {
          CAPTURE(cap);
          const auto expected = computeLocalCostOfClustering(
              enforceMaxClusterId(clusterIds, cap), cellCosts, rate, wiggle, minimalTimestep);
          REQUIRE(batched[cap] == expected);
        }
      }
    }
  }
}

TEST_CASE("Auto merging of clusters") {
  using namespace seissol::initializer::time_stepping;
  const auto clusterIds = std::vector<int>{0, 0, 0, 0, 1, 1, 2};
  const auto cellCosts = std::vector<int>{1, 1, 1, 1, 3, 3, 9};
  const auto minDt = 0.5;
  const auto costBeforeRate2 = computeLocalCostOfClustering(clusterIds, cellCosts, {2}, 1.0, minDt);
  const auto costBeforeRate3 = computeLocalCostOfClustering(clusterIds, cellCosts, {3}, 1.0, minDt);

  SUBCASE("Reduces to GTS") {
    const auto should = 0;
    const auto is = computeMaxClusterIdAfterAutoMerge(
        clusterIds, cellCosts, {2}, std::numeric_limits<double>::max(), 1.0, minDt);
    REQUIRE(is == should);
  }

  SUBCASE("Does nothing for GTS") {
    const auto should = 0;
    for (int i = 1; i <= 5; ++i) {
      const auto is = computeMaxClusterIdAfterAutoMerge(
          enforceMaxClusterId(clusterIds, 0), cellCosts, {1}, i, 0, 0);
      REQUIRE(is == should);
    }
  }

  SUBCASE("No performance loss allowed") {
    const auto should = 2;
    SUBCASE("Rate 2") {
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, {2}, costBeforeRate2, 1, minDt);
      REQUIRE(is == should);
    }
    SUBCASE("Rate 3") {
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, {3}, costBeforeRate3, 1, minDt);
      REQUIRE(is == should);
    }
  }

  SUBCASE("Some performance loss allowed") {
    SUBCASE("Merge one cluster") {
      const auto should = 1;
      const auto is = computeMaxClusterIdAfterAutoMerge(
          clusterIds, cellCosts, {2}, 1.25 * costBeforeRate2, 1, minDt);
      REQUIRE(is == should);
    }
    SUBCASE("Merge two clusters") {
      const auto should = 0;
      const auto is = computeMaxClusterIdAfterAutoMerge(
          clusterIds, cellCosts, {2}, 2.06 * costBeforeRate2, 1, minDt);
      REQUIRE(is == should);
    }
  }
}

TEST_CASE("LTS clustering invariants on a mesh") {
  // Value-free characterization of the end-to-end clustering: rather than pinning golden
  // cluster ids (which would have to be regenerated for every mesh change), this asserts
  // the two structural properties that the refactor must not break.
  std::cout.setstate(std::ios_base::failbit);
  using namespace seissol::initializer::time_stepping;

  // a non-uniform ladder plus an active wiggle sweep -- neither is covered above
  const auto rate = std::vector<std::uint64_t>{4, 2};
  constexpr int MaxClusters = 100;

  const LtsWeightsConfig config{
      seissol::initializer::parameters::BoundaryFormat::I32, rate, 1, 1, 1};

  const seissol::initializer::parameters::LtsParameters ltsParameters(
      rate,
      0.5,   // wiggle factor minimum: enables the sweep
      0.01,  // wiggle factor stepsize
      false, // do not enforce the maximum difference inside the sweep
      MaxClusters,
      false, // no auto merge
      1.0,
      seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor,
      seissol::initializer::parameters::LtsWeightsTypes::ExponentialWeights);

  seissol::initializer::parameters::SeisSolParameters seissolParameters{};
  seissolParameters.timeStepping.lts = ltsParameters;
  seissolParameters.timeStepping.cfl = 1;
  seissolParameters.timeStepping.maxTimestepWidth = 5000.0;
  seissolParameters.model.materialFileName = tpath("Testing/material.yaml");
  seissolParameters.model.useCellHomogenizedMaterial = false;
  seissolParameters.model.plasticity = false;
  const utils::Env env("SEISSOL_");
  seissol::SeisSol seissolInstance(seissolParameters, env);

  auto ltsWeights = std::make_unique<ExponentialWeights>(config, seissolInstance);
  const auto faceMap = defaultFaceMap();
  const auto pumlReader =
      seissol::geometry::PUMLReader(tpath("Testing/mesh.h5"),
                                    "Default",
                                    faceMap,
                                    seissol::initializer::parameters::BoundaryFormat::I32,
                                    seissol::initializer::parameters::TopologyFormat::Geometric,
                                    ltsWeights.get());
  std::cout.clear();

  const auto wiggle = ltsWeights->getWiggleFactor();
  REQUIRE(wiggle > 0.0);
  REQUIRE(wiggle <= 1.0);

  const auto& elements = pumlReader.getElements();
  REQUIRE(!elements.empty());

  double globalMinTimestep = std::numeric_limits<double>::max();
  for (const auto& element : elements) {
    globalMinTimestep = std::min(globalMinTimestep, element.timestep);
  }
  MPI_Allreduce(MPI_IN_PLACE, &globalMinTimestep, 1, MPI_DOUBLE, MPI_MIN, seissol::Mpi::mpi.comm());
  REQUIRE(globalMinTimestep > 0.0);

  // (1) enforceMaximumDifference() and enforceMaxClusterId() only ever move cells to a
  //     finer cluster, so the plain binning is a pointwise upper bound on the result.
  for (const auto& element : elements) {
    CAPTURE(element.globalId);
    const auto unsmoothed = getCluster(element.timestep, globalMinTimestep, wiggle, rate);
    REQUIRE(element.clusterId >= 0);
    REQUIRE(static_cast<std::uint64_t>(element.clusterId) <= unsmoothed);
    REQUIRE(element.clusterId < MaxClusters);
  }

  // (1b) the published ladder describes exactly the clustering that was produced
  const auto& effectiveRates = ltsWeights->effectiveRates();
  for (std::size_t k = 1; k < effectiveRates.size(); ++k) {
    REQUIRE(effectiveRates[k] != 1);
  }
  {
    double globalMaxTimestep = 0.0;
    for (const auto& element : elements) {
      globalMaxTimestep = std::max(globalMaxTimestep, element.timestep);
    }
    MPI_Allreduce(
        MPI_IN_PLACE, &globalMaxTimestep, 1, MPI_DOUBLE, MPI_MAX, seissol::Mpi::mpi.comm());

    const auto published = seissol::initializer::ClusterLadder::forBinning(
        effectiveRates, globalMinTimestep, wiggle, globalMaxTimestep);
    for (const auto& element : elements) {
      CAPTURE(element.globalId);
      REQUIRE(static_cast<std::size_t>(element.clusterId) <= effectiveRates.size());
      // normalizing must not move a single cell
      REQUIRE(published.clusterOf(element.timestep) ==
              getCluster(element.timestep, globalMinTimestep, wiggle, rate));
    }
  }

  // (2) the index-difference-of-one property. This is not a heuristic: the ghost cluster
  //     construction asserts it at runtime, cf. the asserts in TimeManager::addClusters().
  const auto rank = seissol::Mpi::mpi.rank();
  for (const auto& element : elements) {
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto faceType = static_cast<FaceType>(element.boundaries[f]);
      if (!isInternalFaceType(faceType)) {
        continue;
      }
      if (element.neighborRanks[f] != rank) {
        continue; // ghost neighbour, not resolvable from the local element list
      }
      const auto neighbor = static_cast<std::size_t>(element.neighbors[f]);
      if (neighbor >= elements.size()) {
        continue; // domain boundary sentinel
      }
      CAPTURE(element.globalId);
      const auto difference = element.clusterId > elements[neighbor].clusterId
                                  ? element.clusterId - elements[neighbor].clusterId
                                  : elements[neighbor].clusterId - element.clusterId;
      if (faceType == FaceType::DynamicRupture) {
        REQUIRE(difference == 0);
      } else {
        REQUIRE(difference <= 1);
      }
    }
  }
}

} // namespace seissol::unit_test
