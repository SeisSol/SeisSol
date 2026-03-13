// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Reader/Distribution.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

namespace seissol::unit_test {
using seissol::io::reader::Distributor;

// ---------------------------------------------------------------------------
// Distributor: identity mapping (each rank keeps its own data)
// ---------------------------------------------------------------------------

TEST_CASE("Distributor identity mapping" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();

  constexpr std::size_t elemsPerRank = 4;
  const std::size_t offset = rank * elemsPerRank;

  // Each rank owns IDs [offset, offset+1, ..., offset+K-1]
  std::vector<std::size_t> sourceIds(elemsPerRank);
  std::iota(sourceIds.begin(), sourceIds.end(), offset);

  // Each rank requests the same IDs it owns (identity)
  std::vector<std::size_t> targetIds = sourceIds;

  // Source data: value = globalId * 10
  std::vector<double> sourceData(elemsPerRank);
  for (std::size_t i = 0; i < elemsPerRank; ++i) {
    sourceData[i] = static_cast<double>(sourceIds[i]) * 10.0;
  }

  Distributor dist(MPI_COMM_WORLD);
  dist.setup(sourceIds, targetIds);

  std::vector<double> targetData(elemsPerRank, -1.0);
  dist.distribute(targetData.data(), sourceData.data()).complete();

  for (std::size_t i = 0; i < elemsPerRank; ++i) {
    CHECK(targetData[i] == doctest::Approx(targetIds[i] * 10.0));
  }
}

// ---------------------------------------------------------------------------
// Distributor: reversed mapping (each rank gets data from "opposite" rank)
// ---------------------------------------------------------------------------

TEST_CASE("Distributor cross-rank exchange" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  constexpr std::size_t elemsPerRank = 3;

  // Each rank owns IDs [rank*K, rank*K+K-1]
  std::vector<std::size_t> sourceIds(elemsPerRank);
  std::iota(sourceIds.begin(), sourceIds.end(), rank * elemsPerRank);

  // Each rank requests IDs from the "mirror" rank:
  // rank 0 → last rank's IDs, rank 1 → second-to-last, etc.
  const int mirrorRank = size - 1 - rank;
  std::vector<std::size_t> targetIds(elemsPerRank);
  std::iota(targetIds.begin(), targetIds.end(), mirrorRank * elemsPerRank);

  // Source data: value = globalId * 7 + 3
  std::vector<double> sourceData(elemsPerRank);
  for (std::size_t i = 0; i < elemsPerRank; ++i) {
    sourceData[i] = static_cast<double>(sourceIds[i]) * 7.0 + 3.0;
  }

  Distributor dist(MPI_COMM_WORLD);
  dist.setup(sourceIds, targetIds);

  std::vector<double> targetData(elemsPerRank, -1.0);
  dist.distribute(targetData.data(), sourceData.data()).complete();

  // After distribution, targetData[i] should correspond to targetIds[i]
  for (std::size_t i = 0; i < elemsPerRank; ++i) {
    double expected = static_cast<double>(targetIds[i]) * 7.0 + 3.0;
    CHECK(targetData[i] == doctest::Approx(expected));
  }
}

// ---------------------------------------------------------------------------
// Distributor: single element per rank
// ---------------------------------------------------------------------------

TEST_CASE("Distributor single element" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  // Each rank owns 1 element
  std::vector<std::size_t> sourceIds = {static_cast<std::size_t>(rank)};

  // Each rank requests from the next rank (circular)
  std::vector<std::size_t> targetIds = {static_cast<std::size_t>((rank + 1) % size)};

  int sourceData = rank * 1000;

  Distributor dist(MPI_COMM_WORLD);
  dist.setup(sourceIds, targetIds);

  int targetData = -1;
  dist.distribute(&targetData, &sourceData).complete();

  int expectedRank = (rank + 1) % size;
  CHECK(targetData == expectedRank * 1000);
}

} // namespace seissol::unit_test
