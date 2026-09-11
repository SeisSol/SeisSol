// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/InitProcedure/Internal/LtsSetup.h"
#include "Initializer/LtsSetup.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol::unit_test {

namespace {

constexpr std::array<BufferType, BufferCount> AllBufferTypes{
    BufferType::StepIntegrals, BufferType::Derivatives, BufferType::AccumulatedIntegrals};

struct CellSetup {
  std::array<FaceType, Cell::NumFaces> faceTypes{
      FaceType::Regular, FaceType::Regular, FaceType::Regular, FaceType::Regular};
  std::array<std::uint64_t, Cell::NumFaces> neighborClusters{};
  std::uint64_t clusterId{};
};

auto derive(const CellSetup& setup) -> LtsSetup {
  CellLocalInformation primary{};
  primary.faceTypes = setup.faceTypes;
  SecondaryCellLocalInformation secondary{};
  secondary.clusterId = setup.clusterId;
  return initializer::internal::getLtsSetup(primary, secondary, setup.neighborClusters);
}

} // namespace

TEST_CASE("LtsSetup fields do not overlap" * doctest::test_suite("initializer")) {
  LtsSetup setup;

  CHECK_FALSE(setup.hasAnyBuffer());
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    CHECK(setup.neighborBuffer(face) == BufferType::StepIntegrals);
    CHECK_FALSE(setup.neighborGTSRelation(face));
  }

  // fill every field with its largest value; if any two fields shared a bit, at least one of the
  // read-backs below would come out wrong
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    setup.setNeighborBuffer(face, BufferType::AccumulatedIntegrals);
    setup.setNeighborGTSRelation(face, true);
  }
  for (const auto type : AllBufferTypes) {
    setup.setHasBuffer(true, type);
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    CHECK(setup.neighborBuffer(face) == BufferType::AccumulatedIntegrals);
    CHECK(setup.neighborGTSRelation(face));
  }
  for (const auto type : AllBufferTypes) {
    CHECK(setup.hasBuffer(type));
  }
  CHECK(setup.hasAnyBuffer());
}

TEST_CASE("LtsSetup neighbor buffer types are per-face" * doctest::test_suite("initializer")) {
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (const auto type : AllBufferTypes) {
      CAPTURE(face);
      CAPTURE(static_cast<int>(type));

      LtsSetup setup;
      setup.setNeighborBuffer(face, type);

      for (std::size_t other = 0; other < Cell::NumFaces; ++other) {
        const auto expected = other == face ? type : BufferType::StepIntegrals;
        CHECK(setup.neighborBuffer(other) == expected);
      }

      // the neighbor field says nothing about what this cell stores
      CHECK_FALSE(setup.hasAnyBuffer());
    }
  }
}

TEST_CASE("LtsSetup overwrites a neighbor buffer type completely" *
          doctest::test_suite("initializer")) {
  // a set-without-clear would leave the bits of the previous type standing
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    for (const auto type : AllBufferTypes) {
      CAPTURE(face);
      CAPTURE(static_cast<int>(type));

      LtsSetup setup;
      setup.setNeighborBuffer(face, type);
      setup.setNeighborBuffer(face, BufferType::StepIntegrals);
      CHECK(setup.neighborBuffer(face) == BufferType::StepIntegrals);
    }
  }
}

TEST_CASE("LtsSetup storage flags are per-type" * doctest::test_suite("initializer")) {
  for (const auto type : AllBufferTypes) {
    CAPTURE(static_cast<int>(type));

    LtsSetup setup;
    setup.setHasBuffer(true, type);
    for (const auto other : AllBufferTypes) {
      CHECK(setup.hasBuffer(other) == (other == type));
    }
    CHECK(setup.hasAnyBuffer());

    setup.setHasBuffer(false, type);
    CHECK_FALSE(setup.hasBuffer(type));
    CHECK_FALSE(setup.hasAnyBuffer());
  }
}

TEST_CASE("LtsSetup survives the bitmap round trip" * doctest::test_suite("initializer")) {
  // the bitmap is what gets sent to the ghost layers via MPI
  LtsSetup setup;
  setup.setNeighborBuffer(0, BufferType::Derivatives);
  setup.setNeighborBuffer(1, BufferType::AccumulatedIntegrals);
  setup.setNeighborBuffer(2, BufferType::StepIntegrals);
  setup.setNeighborBuffer(3, BufferType::Derivatives);
  setup.setNeighborGTSRelation(1, true);
  setup.setNeighborGTSRelation(2, true);
  setup.setHasBuffer(true, BufferType::StepIntegrals);
  setup.setHasBuffer(true, BufferType::Derivatives);

  const LtsSetup restored(setup.unwrap());

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    CHECK(restored.neighborBuffer(face) == setup.neighborBuffer(face));
    CHECK(restored.neighborGTSRelation(face) == setup.neighborGTSRelation(face));
  }
  for (const auto type : AllBufferTypes) {
    CHECK(restored.hasBuffer(type) == setup.hasBuffer(type));
  }
}

TEST_CASE("LTS setup of a pure GTS cell" * doctest::test_suite("initializer")) {
  CellSetup input;
  input.clusterId = 3;
  input.neighborClusters = {3, 3, 3, 3};

  const auto setup = derive(input);

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    CHECK(setup.neighborBuffer(face) == BufferType::StepIntegrals);
    CHECK(setup.neighborGTSRelation(face));
  }
  CHECK(setup.hasBuffer(BufferType::StepIntegrals));
  CHECK_FALSE(setup.hasBuffer(BufferType::AccumulatedIntegrals));
  CHECK_FALSE(setup.hasBuffer(BufferType::Derivatives));
}

TEST_CASE("LTS setup next to a coarser cluster" * doctest::test_suite("initializer")) {
  CellSetup input;
  input.clusterId = 2;
  input.neighborClusters = {2, 3, 2, 2};

  const auto setup = derive(input);

  // the coarser neighbor has to hand us its full space-time evolution, and we accumulate for it
  CHECK(setup.neighborBuffer(1) == BufferType::Derivatives);
  CHECK_FALSE(setup.neighborGTSRelation(1));
  CHECK(setup.hasBuffer(BufferType::AccumulatedIntegrals));

  for (const auto face : std::array<std::size_t, 3>{0, 2, 3}) {
    CHECK(setup.neighborBuffer(face) == BufferType::StepIntegrals);
    CHECK(setup.neighborGTSRelation(face));
  }
  CHECK(setup.hasBuffer(BufferType::StepIntegrals));

  // accumulating for one neighbor must not force derivatives on everyone else
  CHECK_FALSE(setup.hasBuffer(BufferType::Derivatives));
}

TEST_CASE("LTS setup next to a finer cluster" * doctest::test_suite("initializer")) {
  CellSetup input;
  input.clusterId = 4;
  input.neighborClusters = {3, 4, 4, 4};

  const auto setup = derive(input);

  CHECK(setup.neighborBuffer(0) == BufferType::AccumulatedIntegrals);
  CHECK_FALSE(setup.neighborGTSRelation(0));
  CHECK(setup.hasBuffer(BufferType::Derivatives));
  CHECK(setup.hasBuffer(BufferType::StepIntegrals));
  CHECK_FALSE(setup.hasBuffer(BufferType::AccumulatedIntegrals));
}

TEST_CASE("LTS setup with all three storage types" * doctest::test_suite("initializer")) {
  CellSetup input;
  input.clusterId = 3;
  input.neighborClusters = {2, 3, 4, 3};

  const auto setup = derive(input);

  CHECK(setup.neighborBuffer(0) == BufferType::AccumulatedIntegrals);
  CHECK(setup.neighborBuffer(1) == BufferType::StepIntegrals);
  CHECK(setup.neighborBuffer(2) == BufferType::Derivatives);
  CHECK(setup.neighborBuffer(3) == BufferType::StepIntegrals);

  for (const auto type : AllBufferTypes) {
    CHECK(setup.hasBuffer(type));
  }
}

TEST_CASE("LTS setup ignores external boundary faces" * doctest::test_suite("initializer")) {
  CellSetup input;
  input.clusterId = 1;
  input.neighborClusters = {1, 1, 1, 1};
  input.faceTypes = {
      FaceType::FreeSurface, FaceType::Regular, FaceType::Dirichlet, FaceType::Outflow};

  const auto setup = derive(input);

  CHECK(setup.neighborBuffer(1) == BufferType::StepIntegrals);
  CHECK(setup.neighborGTSRelation(1));

  // boundary faces do not run the neighbor kernel, so they request nothing
  for (const auto face : std::array<std::size_t, 3>{0, 2, 3}) {
    CHECK_FALSE(setup.neighborGTSRelation(face));
  }
  CHECK(setup.hasBuffer(BufferType::StepIntegrals));
  CHECK_FALSE(setup.hasBuffer(BufferType::Derivatives));
  CHECK_FALSE(setup.hasBuffer(BufferType::AccumulatedIntegrals));
}

TEST_CASE("LTS setup of a dynamic rupture face" * doctest::test_suite("initializer")) {
  CellSetup input;
  input.clusterId = 2;
  input.neighborClusters = {2, 2, 2, 2};
  input.faceTypes = {
      FaceType::DynamicRupture, FaceType::Regular, FaceType::Regular, FaceType::Regular};

  const auto setup = derive(input);

  CHECK(setup.neighborBuffer(0) == BufferType::Derivatives);
  CHECK(setup.neighborGTSRelation(0));
  CHECK(setup.hasBuffer(BufferType::Derivatives));
  CHECK(setup.hasBuffer(BufferType::StepIntegrals));
  CHECK_FALSE(setup.hasBuffer(BufferType::AccumulatedIntegrals));
}

TEST_CASE("LTS setup requests only storage that the neighbor keeps" *
          doctest::test_suite("initializer")) {
  // the invariant that is only checked across MPI boundaries today: whatever we request from a
  // face neighbor, that neighbor has to store
  constexpr std::array<std::uint64_t, 3> Clusters{1, 2, 3};

  for (const auto ownCluster : Clusters) {
    for (const auto neighborCluster : Clusters) {
      CAPTURE(ownCluster);
      CAPTURE(neighborCluster);

      CellSetup own;
      own.clusterId = ownCluster;
      own.neighborClusters.fill(neighborCluster);

      CellSetup neighbor;
      neighbor.clusterId = neighborCluster;
      neighbor.neighborClusters.fill(ownCluster);

      const auto ownSetup = derive(own);
      const auto neighborSetup = derive(neighbor);

      CHECK(ownSetup.hasAnyBuffer());
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        CHECK(neighborSetup.hasBuffer(ownSetup.neighborBuffer(face)));
      }
    }
  }
}

TEST_CASE("LTS setup serves a GTS and a coarser neighbor at once" *
          doctest::test_suite("initializer")) {
  // the case that used to need the "GTS on derivatives" normalization: our GTS neighbor also
  // accumulates for a coarser cluster
  CellSetup own;
  own.clusterId = 2;
  own.neighborClusters = {2, 2, 2, 2};

  CellSetup neighbor;
  neighbor.clusterId = 2;
  neighbor.neighborClusters = {2, 3, 2, 2};

  const auto ownSetup = derive(own);
  const auto neighborSetup = derive(neighbor);

  CHECK(ownSetup.neighborBuffer(0) == BufferType::StepIntegrals);
  CHECK(neighborSetup.hasBuffer(ownSetup.neighborBuffer(0)));
  CHECK(neighborSetup.hasBuffer(BufferType::AccumulatedIntegrals));
}

} // namespace seissol::unit_test
