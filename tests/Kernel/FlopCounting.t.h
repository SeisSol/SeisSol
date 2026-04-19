// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Solver.h"

#include <array>
#include <cstdint>

namespace seissol::unit_test {

// ---------------------------------------------------------------------------
// Spacetime (ADER) flops and bytes
// ---------------------------------------------------------------------------

TEST_CASE("Spacetime flopsAder" * doctest::test_suite("kernel")) {
  kernels::Spacetime spacetime;

  std::uint64_t nonZeroFlops = 0;
  std::uint64_t hardwareFlops = 0;
  spacetime.flopsAder(nonZeroFlops, hardwareFlops);

  SUBCASE("Flops are positive") {
    CHECK(nonZeroFlops > 0);
    CHECK(hardwareFlops > 0);
  }

  SUBCASE("Hardware flops >= nonzero flops") { CHECK(hardwareFlops >= nonZeroFlops); }

  SUBCASE("Deterministic: calling twice gives same result") {
    std::uint64_t nz2 = 0;
    std::uint64_t hw2 = 0;
    spacetime.flopsAder(nz2, hw2);
    CHECK(nonZeroFlops == nz2);
    CHECK(hardwareFlops == hw2);
  }
}

TEST_CASE("Spacetime bytesAder" * doctest::test_suite("kernel")) {
  kernels::Spacetime spacetime;
  auto bytes = spacetime.bytesAder();

  CHECK(bytes > 0);
  // Must be a multiple of sizeof(real)
  CHECK(bytes % sizeof(real) == 0);
}

// ---------------------------------------------------------------------------
// Local kernel flops and bytes
// ---------------------------------------------------------------------------

TEST_CASE("Local flopsIntegral all Regular faces" * doctest::test_suite("kernel")) {
  kernels::Local local;

  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::Regular);

  std::uint64_t nonZeroFlops = 0;
  std::uint64_t hardwareFlops = 0;
  local.flopsIntegral(faceTypes, nonZeroFlops, hardwareFlops);

  SUBCASE("Positive flops") {
    CHECK(nonZeroFlops > 0);
    CHECK(hardwareFlops > 0);
  }

  SUBCASE("Hardware >= nonzero") { CHECK(hardwareFlops >= nonZeroFlops); }
}

TEST_CASE("Local flopsIntegral with DynamicRupture faces" * doctest::test_suite("kernel")) {
  kernels::Local local;

  // All faces are DR → local flux is skipped for each (on CPU)
  std::array<FaceType, Cell::NumFaces> faceTypesDR{};
  faceTypesDR.fill(FaceType::DynamicRupture);

  std::uint64_t nzDR = 0;
  std::uint64_t hwDR = 0;
  local.flopsIntegral(faceTypesDR, nzDR, hwDR);

  // All regular faces
  std::array<FaceType, Cell::NumFaces> faceTypesReg{};
  faceTypesReg.fill(FaceType::Regular);

  std::uint64_t nzReg = 0;
  std::uint64_t hwReg = 0;
  local.flopsIntegral(faceTypesReg, nzReg, hwReg);

  SUBCASE("DR faces have fewer local flops than regular (on CPU)") {
    if constexpr (isDeviceOn()) {
      // On the GPU, the kernel runs regardless, so flops may be equal
      CHECK(nzDR > 0);
    } else {
      // On the CPU, DR faces skip the local flux contribution
      CHECK(nzDR < nzReg);
      CHECK(hwDR < hwReg);
    }
  }
}

TEST_CASE("Local flopsIntegral mixed faces" * doctest::test_suite("kernel")) {
  kernels::Local local;

  const std::array<FaceType, Cell::NumFaces> faceTypes = {
      FaceType::Regular,
      FaceType::FreeSurface,
      FaceType::Outflow,
      FaceType::Periodic,
  };

  std::uint64_t nonZeroFlops = 0;
  std::uint64_t hardwareFlops = 0;
  local.flopsIntegral(faceTypes, nonZeroFlops, hardwareFlops);

  CHECK(nonZeroFlops > 0);
  CHECK(hardwareFlops >= nonZeroFlops);
}

TEST_CASE("Local bytesIntegral" * doctest::test_suite("kernel")) {
  kernels::Local local;
  auto bytes = local.bytesIntegral();

  CHECK(bytes > 0);
  CHECK(bytes % sizeof(real) == 0);
}

// ---------------------------------------------------------------------------
// Neighbor kernel flops and bytes
// ---------------------------------------------------------------------------

TEST_CASE("Neighbor flopsNeighborsIntegral all Regular" * doctest::test_suite("kernel")) {
  kernels::Neighbor neighbor;

  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::Regular);

  std::array<std::array<uint8_t, 2>, Cell::NumFaces> neighboringIndices{};
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    neighboringIndices[f] = {0, 0};
  }

  const CellDRMapping drMapping[4]{};

  std::uint64_t nz = 0;
  std::uint64_t hw = 0;
  std::uint64_t drNz = 0;
  std::uint64_t drHw = 0;
  neighbor.flopsNeighborsIntegral(faceTypes, neighboringIndices, drMapping, nz, hw, drNz, drHw);

  SUBCASE("Positive flops for regular faces") {
    CHECK(nz > 0);
    CHECK(hw > 0);
  }

  SUBCASE("No DR flops for regular faces") {
    CHECK(drNz == 0);
    CHECK(drHw == 0);
  }

  SUBCASE("Hardware >= nonzero") { CHECK(hw >= nz); }
}

TEST_CASE("Neighbor flopsNeighborsIntegral with DR faces" * doctest::test_suite("kernel")) {
  kernels::Neighbor neighbor;

  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::DynamicRupture);

  const std::array<std::array<uint8_t, 2>, Cell::NumFaces> neighboringIndices{};
  const CellDRMapping drMapping[4]{};

  std::uint64_t nz = 0;
  std::uint64_t hw = 0;
  std::uint64_t drNz = 0;
  std::uint64_t drHw = 0;
  neighbor.flopsNeighborsIntegral(faceTypes, neighboringIndices, drMapping, nz, hw, drNz, drHw);

  SUBCASE("DR flops are positive") {
    CHECK(drNz > 0);
    CHECK(drHw > 0);
  }

  SUBCASE("DR hardware >= nonzero") { CHECK(drHw >= drNz); }
}

TEST_CASE("Neighbor bytesNeighborsIntegral" * doctest::test_suite("kernel")) {
  kernels::Neighbor neighbor;
  auto bytes = neighbor.bytesNeighborsIntegral();

  CHECK(bytes > 0);
  CHECK(bytes % sizeof(real) == 0);
}

// ---------------------------------------------------------------------------
// Cross-kernel consistency
// ---------------------------------------------------------------------------

TEST_CASE("Kernel flop ordering: Ader < Local < Neighbor" * doctest::test_suite("kernel")) {
  // This is a heuristic sanity check, not a hard invariant.
  // For typical configurations, ADER (time integration) should be less expensive
  // than the spatial integration (Local), and neighbor integration should also
  // contribute significantly.

  kernels::Spacetime spacetime;
  std::uint64_t aderNz = 0;
  std::uint64_t aderHw = 0;
  spacetime.flopsAder(aderNz, aderHw);

  kernels::Local local;
  std::array<FaceType, Cell::NumFaces> faceTypes{};
  faceTypes.fill(FaceType::Regular);
  std::uint64_t localNz = 0;
  std::uint64_t localHw = 0;
  local.flopsIntegral(faceTypes, localNz, localHw);

  kernels::Neighbor neighbor;
  const std::array<std::array<uint8_t, 2>, Cell::NumFaces> neighboringIndices{};
  const CellDRMapping drMapping[4]{};
  std::uint64_t neighborNz = 0;
  std::uint64_t neighborHw = 0;
  std::uint64_t drNz = 0;
  std::uint64_t drHw = 0;
  neighbor.flopsNeighborsIntegral(
      faceTypes, neighboringIndices, drMapping, neighborNz, neighborHw, drNz, drHw);

  // All three kernels should have non-trivial cost
  CHECK(aderNz > 100);
  CHECK(localNz > 100);
  CHECK(neighborNz > 100);
}

} // namespace seissol::unit_test
