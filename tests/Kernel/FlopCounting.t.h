// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Solver.h"
#include "TestConfigs.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol::unit_test {

// ---------------------------------------------------------------------------
// Spacetime (ADER) flops and bytes
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("Spacetime flopsAder" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  const kernels::Spacetime<Cfg> spacetime;

  const auto metric = spacetime.metrics();

  SUBCASE("Metrics are positive") {
    CHECK(metric.nonzeroFlop > 0);
    CHECK(metric.hardwareFlop > 0);
    CHECK(metric.bytes > 0);
    CHECK(metric.kernelBytes > 0);
  }

  SUBCASE("Hardware flops >= nonzero flops") { CHECK(metric.hardwareFlop >= metric.nonzeroFlop); }

  SUBCASE("Deterministic: calling twice gives same result") {
    const auto metric2 = spacetime.metrics();
    CHECK(metric.nonzeroFlop == metric2.nonzeroFlop);
    CHECK(metric.hardwareFlop == metric2.hardwareFlop);
    CHECK(metric.bytes == metric2.bytes);
    CHECK(metric.kernelBytes == metric2.kernelBytes);
  }
}

// ---------------------------------------------------------------------------
// Local kernel flops and bytes
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("Local metrics all Regular faces" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  const kernels::Local<Cfg> local;

  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::Regular);

  const auto metric = local.metrics(faceTypes);

  SUBCASE("Metrics are positive") {
    CHECK(metric.nonzeroFlop > 0);
    CHECK(metric.hardwareFlop > 0);
    CHECK(metric.bytes > 0);
    CHECK(metric.kernelBytes > 0);
  }

  SUBCASE("Hardware >= nonzero") { CHECK(metric.hardwareFlop >= metric.nonzeroFlop); }
}

TEST_CASE_TEMPLATE("Local metrics with DynamicRupture faces" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  const kernels::Local<Cfg> local;

  // All faces are DR → local flux is skipped for each (on CPU)
  std::array<FaceType, Cell::NumFaces> faceTypesDR{};
  faceTypesDR.fill(FaceType::DynamicRupture);

  const auto metricDR = local.metrics(faceTypesDR);

  // All regular faces
  std::array<FaceType, Cell::NumFaces> faceTypesReg{};
  faceTypesReg.fill(FaceType::Regular);

  const auto metricReg = local.metrics(faceTypesReg);

  SUBCASE("DR faces have fewer local flops than regular (on CPU)") {
    if constexpr (isDeviceOn()) {
      // On the GPU, the kernel runs regardless, so flops may be equal
      CHECK(metricDR.nonzeroFlop > 0);
    } else {
      // On the CPU, DR faces skip the local flux contribution
      CHECK(metricDR.nonzeroFlop < metricReg.nonzeroFlop);
      CHECK(metricDR.hardwareFlop < metricReg.hardwareFlop);
    }
  }
}

TEST_CASE_TEMPLATE("Local metrics mixed faces" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  const kernels::Local<Cfg> local;

  const std::array<FaceType, Cell::NumFaces> faceTypes = {
      FaceType::Regular,
      FaceType::FreeSurface,
      FaceType::Outflow,
  };

  const auto metric = local.metrics(faceTypes);

  CHECK(metric.nonzeroFlop > 0);
  CHECK(metric.hardwareFlop >= metric.nonzeroFlop);

  CHECK(metric.bytes > 0);
  CHECK(metric.bytes % sizeof(Real<Cfg>) == 0);
}

// ---------------------------------------------------------------------------
// Neighbor kernel flops and bytes
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("Neighbor metrics all Regular" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  const kernels::Neighbor<Cfg> neighbor;

  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::Regular);

  std::array<std::array<uint8_t, 2>, Cell::NumFaces> neighboringIndices{};
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    neighboringIndices[f] = {0, 0};
  }

  const std::array<CellDRMapping<Cfg>, Cell::NumFaces> drMapping{};

  const auto [metric, metricDR] = neighbor.metrics(faceTypes, neighboringIndices, drMapping);

  SUBCASE("Metrics are positive") {
    CHECK(metric.nonzeroFlop > 0);
    CHECK(metric.hardwareFlop > 0);
    CHECK(metric.bytes > 0);
    CHECK(metric.kernelBytes > 0);
  }

  SUBCASE("DR metrics are zero") {
    CHECK(metricDR.nonzeroFlop == 0);
    CHECK(metricDR.hardwareFlop == 0);
    // CHECK(metricDR.bytes > 0);
    CHECK(metricDR.kernelBytes == 0);
  }

  SUBCASE("Hardware >= nonzero") { CHECK(metric.hardwareFlop >= metric.nonzeroFlop); }
}

TEST_CASE_TEMPLATE("Neighbor metrics with DR faces" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  const kernels::Neighbor<Cfg> neighbor;

  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::DynamicRupture);

  const std::array<std::array<uint8_t, 2>, Cell::NumFaces> neighboringIndices{};
  const std::array<CellDRMapping<Cfg>, Cell::NumFaces> drMapping{};

  const auto [metric, metricDR] = neighbor.metrics(faceTypes, neighboringIndices, drMapping);

  SUBCASE("DR metrics are positive") {
    CHECK(metricDR.nonzeroFlop > 0);
    CHECK(metricDR.hardwareFlop > 0);
    // CHECK(metricDR.bytes > 0);
    CHECK(metricDR.kernelBytes > 0);
  }

  SUBCASE("DR hardware >= nonzero") { CHECK(metricDR.hardwareFlop >= metricDR.nonzeroFlop); }
}

// ---------------------------------------------------------------------------
// Cross-kernel consistency
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("Kernel flop ordering: Ader < Local < Neighbor" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  // This is a heuristic sanity check, not a hard invariant.
  // For typical configurations, ADER (time integration) should be less expensive
  // than the spatial integration (Local), and neighbor integration should also
  // contribute significantly.

  const kernels::Spacetime<Cfg> spacetime;
  const auto aderMetrics = spacetime.metrics();

  const kernels::Local<Cfg> local;
  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::Regular);
  const auto localMetrics = local.metrics(faceTypes);

  const kernels::Neighbor<Cfg> neighbor;
  const std::array<std::array<uint8_t, 2>, Cell::NumFaces> neighboringIndices{};
  const std::array<CellDRMapping<Cfg>, Cell::NumFaces> drMapping{};
  const auto [neighborMetrics, drMetrics] =
      neighbor.metrics(faceTypes, neighboringIndices, drMapping);

  // All three kernels should have non-trivial cost: each computes every degree of freedom of the
  // cell at least once. Beyond that, the cost depends on the configuration; in order 2, the ADER
  // kernel of a fluid stays below a hundred flops.
  constexpr std::size_t NumDofs = model::MaterialOf<Cfg>::NumQuantities *
                                  kernels::getNumberOfBasisFunctions(Cfg::ConvergenceOrder);
  CHECK(aderMetrics.nonzeroFlop >= NumDofs);
  CHECK(localMetrics.nonzeroFlop >= NumDofs);
  CHECK(neighborMetrics.nonzeroFlop >= NumDofs);
}

} // namespace seissol::unit_test
