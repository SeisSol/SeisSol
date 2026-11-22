// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Constants.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Geometry/FaceMap.h"
#include "Core/Geometry/MeshParameters.h"
#include "Core/Geometry/PUMLReader.h"
#include "Core/IO/Easi/ReaderBuilder.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Parallel/MPI.h"
#include "Core/TimeStepping/Clustering/ClusterLadder.h"
#include "Core/TimeStepping/Clustering/Clustering.h"
#include "Core/TimeStepping/Clustering/ClusteringCost.h"
#include "Core/TimeStepping/Clustering/VertexWeights/WeightsFactory.h"
#include "Core/TimeStepping/Clustering/VertexWeights/WeightsModels.h"
#include "Core/TimeStepping/Layout/GlobalTimestep.h"
#include "Core/TimeStepping/LtsParameters.h"
#include "Simulation/Parameters/SeisSolParameters.h"
#include "TestHelper.h"
#include "WavePropagation/Setup/Timesteps.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::unit_test {

/// The time steps of `cells` with the materials of `parameters`, as the run computes them.
inline seissol::initializer::GlobalTimestep
    timestepsOf(const seissol::initializer::CellToVertexArray& cells,
                const seissol::initializer::parameters::SeisSolParameters& parameters) {
  return seissol::initializer::computeTimesteps(cells,
                                                parameters.runConfigs,
                                                parameters.wavePropagation,
                                                parameters.timeStepping,
                                                parameters.model.materialFileNames,
                                                data::modelLoader());
}

/// The time steps of the cells with the materials of `parameters`, as the run computes them.
inline seissol::initializer::TimestepFunction
    materialTimesteps(const seissol::initializer::parameters::SeisSolParameters& parameters) {
  return [&parameters](const seissol::initializer::CellToVertexArray& cells) {
    return timestepsOf(cells, parameters);
  };
}

TEST_CASE("LTS Weights" * doctest::test_suite("initializer")) {
  std::cout.setstate(std::ios_base::failbit);
  using namespace seissol::initializer;

  const auto faceMap = defaultFaceMap();

  const ClusteringConfig config{.rate = {2},
                                .vertexWeightElement = 1,
                                .vertexWeightDynamicRupture = 1,
                                .vertexWeightFreeSurfaceWithGravity = 1,
                                .faceMap = &faceMap,
                                .configCostFactors = {}};

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
  seissolParameters.model.materialFileNames = {tpath("Testing/material.yaml")};
  seissolParameters.wavePropagation.useCellHomogenizedMaterial = false;
  seissolParameters.plasticity.enabled = false;

  Clustering clustering(config,
                        seissolParameters.runConfigs,
                        seissolParameters.timeStepping.lts,
                        materialTimesteps(seissolParameters));
  ExponentialWeights weightModel;
  const auto pumlReader =
      seissol::geometry::PUMLReader(tpath("Testing/mesh.h5"),
                                    "Default",
                                    faceMap,
                                    seissol::initializer::parameters::BoundaryFormat::I32,
                                    seissol::initializer::parameters::TopologyFormat::Geometric,
                                    &clustering,
                                    &weightModel);
  std::cout.clear();

  const auto givenWeights =
      std::vector<std::uint64_t>(weightModel.vertexWeights(), weightModel.vertexWeights() + 24);

  const auto expectedWeights = std::vector<std::uint64_t>{2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2,
                                                          2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1};

  CHECK(givenWeights == expectedWeights);
}

TEST_CASE("The cost of a cell follows its configuration" * doctest::test_suite("initializer")) {
  // mesh.h5 has its 24 cells in group 1, and no faces with a weight of their own
  using namespace seissol::initializer;

  const auto faceMap = defaultFaceMap();

  const auto cellCosts = [&](std::vector<double> configCostFactors,
                             std::unordered_map<int, ConfigId> groupConfigs) {
    std::cout.setstate(std::ios_base::failbit);
    const ClusteringConfig config{.rate = {2},
                                  .vertexWeightElement = 10,
                                  .vertexWeightDynamicRupture = 1,
                                  .vertexWeightFreeSurfaceWithGravity = 1,
                                  .faceMap = &faceMap,
                                  .configCostFactors = std::move(configCostFactors)};

    seissol::initializer::parameters::SeisSolParameters seissolParameters{};
    seissolParameters.timeStepping.lts = seissol::initializer::parameters::LtsParameters(
        {2},
        1.0,
        0.01,
        false,
        100,
        false,
        1.0,
        seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor,
        seissol::initializer::parameters::LtsWeightsTypes::ExponentialWeights);
    seissolParameters.timeStepping.cfl = 1;
    seissolParameters.timeStepping.maxTimestepWidth = 5000.0;
    seissolParameters.model.materialFileNames = {tpath("Testing/material.yaml")};
    seissolParameters.wavePropagation.useCellHomogenizedMaterial = false;
    seissolParameters.plasticity.enabled = false;
    seissolParameters.runConfigs.groupConfigs = std::move(groupConfigs);

    Clustering clustering(config,
                          seissolParameters.runConfigs,
                          seissolParameters.timeStepping.lts,
                          materialTimesteps(seissolParameters));
    ExponentialWeights weightModel;
    const auto pumlReader =
        seissol::geometry::PUMLReader(tpath("Testing/mesh.h5"),
                                      "Default",
                                      faceMap,
                                      seissol::initializer::parameters::BoundaryFormat::I32,
                                      seissol::initializer::parameters::TopologyFormat::Geometric,
                                      &clustering,
                                      &weightModel);
    std::cout.clear();
    return clustering.result().cellCosts;
  };
  const auto allCells = [](std::uint64_t cost) { return std::vector<std::uint64_t>(24, cost); };

  // without factors, every cell has the element weight
  CHECK(cellCosts({}, {}) == allCells(10));

  // the factor of its configuration scales the element weight of a cell, rounded, but to at least 1
  std::vector<double> factors(builtConfigCount(), 1.0);
  factors[0] = 2.5;
  CHECK(cellCosts(factors, {}) == allCells(25));
  factors[0] = 1.26;
  CHECK(cellCosts(factors, {}) == allCells(13));
  factors[0] = 0.01;
  CHECK(cellCosts(factors, {}) == allCells(1));

  if (builtConfigCount() > 1) {
    // a cell has the configuration of its group
    factors[0] = 1.0;
    factors[1] = 3.0;
    CHECK(cellCosts(factors, {}) == allCells(10));
    CHECK(cellCosts(factors, {{1, 1}}) == allCells(30));
    CHECK(cellCosts(factors, {{2, 1}}) == allCells(10));
  }
}

TEST_CASE("The LTS clustering sees the vertex order of the simulation" *
          doctest::test_suite("initializer")) {
  // The clustering runs before PUMLReader::getMesh(), but has to compute the time steps with the
  // canonical vertex order that getMesh() hands on to the simulation, not with the order of the
  // mesh file: the time step (via round-off) and the cell-homogenized material (via the
  // placement of the quadrature points) both depend on it. mesh.h5 lists most cells in a
  // non-canonical order, and the material varies inside the cells. (Materials without an
  // averaging rule, such as the anisotropic and the poroelastic ones, sample the barycenter
  // either way.)
  using namespace seissol::initializer;

  const auto faceMap = defaultFaceMap();
  const ClusteringConfig config{.rate = {2},
                                .vertexWeightElement = 1,
                                .vertexWeightDynamicRupture = 1,
                                .vertexWeightFreeSurfaceWithGravity = 1,
                                .faceMap = &faceMap,
                                .configCostFactors = {}};

  for (const bool homogenized : {false, true}) {
    CAPTURE(homogenized);
    std::cout.setstate(std::ios_base::failbit);

    seissol::initializer::parameters::SeisSolParameters seissolParameters{};
    seissolParameters.timeStepping.lts = seissol::initializer::parameters::LtsParameters(
        {2},
        1.0,
        0.01,
        false,
        100,
        false,
        1.0,
        seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor,
        seissol::initializer::parameters::LtsWeightsTypes::ExponentialWeights);
    seissolParameters.timeStepping.cfl = 1;
    seissolParameters.timeStepping.maxTimestepWidth = 5000.0;
    seissolParameters.model.materialFileNames = {tpath("Testing/material-graded.yaml")};
    seissolParameters.wavePropagation.useCellHomogenizedMaterial = homogenized;
    seissolParameters.plasticity.enabled = false;

    Clustering clustering(config,
                          seissolParameters.runConfigs,
                          seissolParameters.timeStepping.lts,
                          materialTimesteps(seissolParameters));
    ExponentialWeights weightModel;
    const auto pumlReader =
        seissol::geometry::PUMLReader(tpath("Testing/mesh.h5"),
                                      "Default",
                                      faceMap,
                                      seissol::initializer::parameters::BoundaryFormat::I32,
                                      seissol::initializer::parameters::TopologyFormat::Geometric,
                                      &clustering,
                                      &weightModel);
    const auto recomputed =
        timestepsOf(CellToVertexArray::fromMeshReader(pumlReader), seissolParameters);
    std::cout.clear();

    const auto& elements = pumlReader.getElements();
    REQUIRE(elements.size() == recomputed.cellTimeStepWidths.size());
    for (std::size_t i = 0; i < elements.size(); ++i) {
      CAPTURE(i);
      // bit for bit: the same vertices in the same order give the same floating-point result
      CHECK(elements[i].timestep == recomputed.cellTimeStepWidths[i]);
    }
  }
}

TEST_CASE("A cell takes the smallest time step of the materials of the simulations" *
          doctest::test_suite("initializer")) {
  // The copies of a cell, which hold the simulations of a run, step together, so that the time
  // step of a cell has to suit the materials of all simulations of the run.
  using namespace seissol::initializer;

  const auto faceMap = defaultFaceMap();
  const ClusteringConfig config{.rate = {2},
                                .vertexWeightElement = 1,
                                .vertexWeightDynamicRupture = 1,
                                .vertexWeightFreeSurfaceWithGravity = 1,
                                .faceMap = &faceMap,
                                .configCostFactors = {}};

  const auto timesteps = [&](const std::vector<std::string>& materialFileNames) {
    seissol::initializer::parameters::SeisSolParameters seissolParameters{};
    seissolParameters.timeStepping.lts = seissol::initializer::parameters::LtsParameters(
        {2},
        1.0,
        0.01,
        false,
        100,
        false,
        1.0,
        seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor,
        seissol::initializer::parameters::LtsWeightsTypes::ExponentialWeights);
    seissolParameters.timeStepping.cfl = 1;
    seissolParameters.timeStepping.maxTimestepWidth = 5000.0;
    seissolParameters.model.materialFileNames = materialFileNames;
    seissolParameters.runConfigs.simulationCount = materialFileNames.size();
    seissolParameters.plasticity.enabled = false;

    Clustering clustering(config,
                          seissolParameters.runConfigs,
                          seissolParameters.timeStepping.lts,
                          materialTimesteps(seissolParameters));
    ExponentialWeights weightModel;
    const auto pumlReader =
        seissol::geometry::PUMLReader(tpath("Testing/mesh.h5"),
                                      "Default",
                                      faceMap,
                                      seissol::initializer::parameters::BoundaryFormat::I32,
                                      seissol::initializer::parameters::TopologyFormat::Geometric,
                                      &clustering,
                                      &weightModel);
    return timestepsOf(CellToVertexArray::fromMeshReader(pumlReader), seissolParameters)
        .cellTimeStepWidths;
  };

  std::cout.setstate(std::ios_base::failbit);
  const auto layered = timesteps({tpath("Testing/material.yaml")});
  const auto graded = timesteps({tpath("Testing/material-graded.yaml")});
  const auto both = timesteps({tpath("Testing/material.yaml"),
                               tpath("Testing/material-graded.yaml"),
                               tpath("Testing/material.yaml")});
  std::cout.clear();

  REQUIRE(graded.size() == layered.size());
  REQUIRE(both.size() == layered.size());
  bool materialsDiffer = false;
  for (std::size_t i = 0; i < layered.size(); ++i) {
    CAPTURE(i);
    CHECK(both[i] == std::min(layered[i], graded[i]));
    materialsDiffer = materialsDiffer || layered[i] != graded[i];
  }
  // otherwise, any of them would do
  CHECK(materialsDiffer);
}

TEST_CASE("Cost function for LTS" * doctest::test_suite("initializer")) {
  const auto eps = 10e-12;
  using namespace seissol::initializer;

  SUBCASE("No clusters") {
    const std::vector<std::size_t> clusterIds = {};
    const std::vector<std::uint64_t> cellCosts = {};
    const auto is = computeLocalCostOfClustering(clusterIds, cellCosts, {2}, 1.0, 1.0);
    const auto should = 0.0;
    CHECK(AbsApprox(is).epsilon(eps) == should);
  }

  SUBCASE("One cluster") {
    const std::vector<std::size_t> clusterIds = {0, 0, 0, 0, 0};
    std::vector<std::uint64_t> cellCosts = {1, 2, 3, 4, 5};
    for (int i = 1; i <= 10; ++i) {
      const auto dt = 1.0 / i;
      for (int j = 1; j <= 10; ++j) {
        const auto wiggleFactor = 1.0 / j;

        const auto is = computeLocalCostOfClustering(clusterIds, cellCosts, {2}, wiggleFactor, dt);

        const auto totalCost = static_cast<double>(
            std::accumulate(cellCosts.begin(), cellCosts.end(), std::uint64_t{0}));
        const auto effectiveDt = dt * wiggleFactor;

        const auto should = totalCost * (1.0 / effectiveDt);
        CHECK(AbsApprox(is).epsilon(eps) == should);
      }
    }
  }

  SUBCASE("Two clusters") {
    const std::vector<std::size_t> clusterIds = {1, 0, 1, 1};
    const std::vector<std::uint64_t> cellCosts = {2, 1, 3, 1};
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
          CHECK(AbsApprox(is).epsilon(eps) == should);
        }
      }
    }
  }

  SUBCASE("Three clusters") {
    const std::vector<std::size_t> clusterIds = {2, 0, 1, 1, 1};
    const std::vector<std::uint64_t> cellCosts = {2, 1, 3, 1, 2};
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
          CHECK(AbsApprox(is).epsilon(eps) == should);
        }
      }
    }
  }
}

TEST_CASE("Enforce max cluster id" * doctest::test_suite("initializer")) {
  using namespace seissol::initializer;
  const auto clusterIds = std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 6, 5, 4, 3, 2, 1, 0};
  SUBCASE("No change") {
    const auto& should = clusterIds;
    const auto is = enforceMaxClusterId(clusterIds, 6);
    CHECK(is == should);
  }
  SUBCASE("Only one cluster") {
    const auto should = std::vector<std::size_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const auto is = enforceMaxClusterId(clusterIds, 0);
    CHECK(is == should);
  }
  SUBCASE("Three clusters") {
    const auto should = std::vector<std::size_t>{0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0};
    const auto is = enforceMaxClusterId(clusterIds, 2);
    CHECK(is == should);
  }
}

TEST_CASE("Batched costs of capped clusterings" * doctest::test_suite("initializer")) {
  using namespace seissol::initializer;
  const auto clusterIds = std::vector<std::size_t>{0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 4};
  const auto cellCosts = std::vector<std::uint64_t>{7, 3, 5, 1, 3, 3, 9, 2, 4, 6, 8, 5};
  const std::size_t maxClusterId = 4;

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
        REQUIRE(batched.size() == maxClusterId + 1);
        for (std::size_t cap = 0; cap <= maxClusterId; ++cap) {
          CAPTURE(cap);
          const auto expected = computeLocalCostOfClustering(
              enforceMaxClusterId(clusterIds, cap), cellCosts, rate, wiggle, minimalTimestep);
          CHECK(batched[cap] == AbsApprox(expected).delta(1e-12));
        }
      }
    }
  }
}

TEST_CASE("Auto merging of clusters" * doctest::test_suite("initializer")) {
  using namespace seissol::initializer;
  const auto clusterIds = std::vector<std::size_t>{0, 0, 0, 0, 1, 1, 2};
  const auto cellCosts = std::vector<std::uint64_t>{1, 1, 1, 1, 3, 3, 9};
  const auto minDt = 0.5;
  const auto costBeforeRate2 = computeLocalCostOfClustering(clusterIds, cellCosts, {2}, 1.0, minDt);
  const auto costBeforeRate3 = computeLocalCostOfClustering(clusterIds, cellCosts, {3}, 1.0, minDt);

  SUBCASE("Reduces to GTS") {
    const auto should = 0;
    const auto is = computeMaxClusterIdAfterAutoMerge(
        clusterIds, cellCosts, {2}, std::numeric_limits<double>::max(), 1.0, minDt);
    CHECK(is == should);
  }

  SUBCASE("Does nothing for GTS") {
    const auto should = 0;
    for (int i = 1; i <= 5; ++i) {
      const auto is = computeMaxClusterIdAfterAutoMerge(
          enforceMaxClusterId(clusterIds, 0), cellCosts, {1}, i, 0, 0);
      CHECK(is == should);
    }
  }

  SUBCASE("No performance loss allowed") {
    const auto should = 2;
    SUBCASE("Rate 2") {
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, {2}, costBeforeRate2, 1, minDt);
      CHECK(is == should);
    }
    SUBCASE("Rate 3") {
      const auto is =
          computeMaxClusterIdAfterAutoMerge(clusterIds, cellCosts, {3}, costBeforeRate3, 1, minDt);
      CHECK(is == should);
    }
  }

  SUBCASE("Some performance loss allowed") {
    SUBCASE("Merge one cluster") {
      const auto should = 1;
      const auto is = computeMaxClusterIdAfterAutoMerge(
          clusterIds, cellCosts, {2}, 1.25 * costBeforeRate2, 1, minDt);
      CHECK(is == should);
    }
    SUBCASE("Merge two clusters") {
      const auto should = 0;
      const auto is = computeMaxClusterIdAfterAutoMerge(
          clusterIds, cellCosts, {2}, 2.06 * costBeforeRate2, 1, minDt);
      CHECK(is == should);
    }
  }
}

TEST_CASE("LTS clustering invariants on a mesh" * doctest::test_suite("initializer")) {
  // Value-free characterization of the end-to-end clustering: rather than pinning golden
  // cluster ids (which would have to be regenerated for every mesh change), this asserts
  // the two structural properties that the refactor must not break.
  std::cout.setstate(std::ios_base::failbit);
  using namespace seissol::initializer;

  // a non-uniform ladder plus an active wiggle sweep -- neither is covered above
  const auto rate = std::vector<std::uint64_t>{4, 2};
  constexpr int MaxClusters = 100;

  const auto faceMap = defaultFaceMap();

  const ClusteringConfig config{.rate = rate,
                                .vertexWeightElement = 1,
                                .vertexWeightDynamicRupture = 1,
                                .vertexWeightFreeSurfaceWithGravity = 1,
                                .faceMap = &faceMap,
                                .configCostFactors = {}};

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
  seissolParameters.model.materialFileNames = {tpath("Testing/material.yaml")};
  seissolParameters.wavePropagation.useCellHomogenizedMaterial = false;
  seissolParameters.plasticity.enabled = false;

  Clustering clustering(config,
                        seissolParameters.runConfigs,
                        seissolParameters.timeStepping.lts,
                        materialTimesteps(seissolParameters));
  ExponentialWeights weightModel;
  const auto pumlReader =
      seissol::geometry::PUMLReader(tpath("Testing/mesh.h5"),
                                    "Default",
                                    faceMap,
                                    seissol::initializer::parameters::BoundaryFormat::I32,
                                    seissol::initializer::parameters::TopologyFormat::Geometric,
                                    &clustering,
                                    &weightModel);
  std::cout.clear();

  const auto wiggle = clustering.result().wiggleFactor;
  CHECK(wiggle > 0.0);
  CHECK(wiggle <= 1.0);

  const auto& elements = pumlReader.getElements();
  CHECK(!elements.empty());

  double globalMinTimestep = std::numeric_limits<double>::max();
  for (const auto& element : elements) {
    globalMinTimestep = std::min(globalMinTimestep, element.timestep);
  }
  globalMinTimestep = seissol::Mpi::mpi.allreduce(globalMinTimestep, MPI_MIN);
  CHECK(globalMinTimestep > 0.0);

  // (1) enforceMaximumDifference() and enforceMaxClusterId() only ever move cells to a
  //     finer cluster, so the plain binning is a pointwise upper bound on the result.
  for (const auto& element : elements) {
    CAPTURE(element.globalId);
    const auto unsmoothed = getCluster(element.timestep, globalMinTimestep, wiggle, rate);
    CHECK(element.clusterId >= 0);
    CHECK(static_cast<std::uint64_t>(element.clusterId) <= unsmoothed);
    CHECK(element.clusterId < MaxClusters);
  }

  // (1b) the published ladder describes exactly the clustering that was produced
  const auto& effectiveRates = clustering.result().ratios;
  for (std::size_t k = 1; k < effectiveRates.size(); ++k) {
    CHECK(effectiveRates[k] != 1);
  }
  {
    // `exact`, not `forBinning`: the published vector is a complete ladder, and re-applying the
    // abbreviated convention to it would append clusters that the clustering never had
    const auto published =
        seissol::initializer::ClusterLadder::exact(effectiveRates, globalMinTimestep, wiggle);
    for (const auto& element : elements) {
      CAPTURE(element.globalId);
      CHECK(element.clusterId <= effectiveRates.size());
      // expanding must not move a single cell
      CHECK(published.clusterOf(element.timestep) ==
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
        continue; // ghost neighbor, not resolvable from the local element list
      }
      if (!element.neighbors[f].hasValue()) {
        continue; // domain boundary
      }
      const auto neighbor = element.neighbors[f].value();
      CAPTURE(element.globalId);
      const auto difference = element.clusterId > elements[neighbor].clusterId
                                  ? element.clusterId - elements[neighbor].clusterId
                                  : elements[neighbor].clusterId - element.clusterId;
      if (faceType == FaceType::DynamicRupture) {
        CHECK(difference == 0);
      } else {
        CHECK(difference <= 1);
      }
    }
  }
}

} // namespace seissol::unit_test
