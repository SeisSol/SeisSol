// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/MPI.h"

#include <cmath>
#include <string>
#include <vector>

namespace seissol::unit_test {

// ---------------------------------------------------------------------------
// collect<T>: MPI_Gather to rank 0
// ---------------------------------------------------------------------------

TEST_CASE("MPI collect double" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  // Each rank sends its rank as a double
  const double value = static_cast<double>(rank) * 10.0;
  auto gathered = mpi.collect(value);

  if (rank == 0) {
    REQUIRE(gathered.size() == static_cast<std::size_t>(size));
    for (int r = 0; r < size; ++r) {
      CHECK(gathered[r] == doctest::Approx(r * 10.0));
    }
  }
}

TEST_CASE("MPI collect int" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  const int value = rank + 100;
  auto gathered = mpi.collect(value);

  if (rank == 0) {
    REQUIRE(gathered.size() == static_cast<std::size_t>(size));
    for (int r = 0; r < size; ++r) {
      CHECK(gathered[r] == r + 100);
    }
  }
}

// ---------------------------------------------------------------------------
// collectContainer: MPI_Gatherv with variable-length containers
// ---------------------------------------------------------------------------

TEST_CASE("MPI collectContainer same size" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  // Each rank sends a vector of 3 doubles: [rank, rank+1, rank+2]
  const std::vector<double> local = {rank * 1.0, rank + 1.0, rank + 2.0};
  auto collected = mpi.collectContainer(local);

  if (rank == 0) {
    REQUIRE(collected.size() == static_cast<std::size_t>(size));
    for (int r = 0; r < size; ++r) {
      REQUIRE(collected[r].size() == 3);
      CHECK(collected[r][0] == doctest::Approx(r * 1.0));
      CHECK(collected[r][1] == doctest::Approx(r + 1.0));
      CHECK(collected[r][2] == doctest::Approx(r + 2.0));
    }
  }
}

TEST_CASE("MPI collectContainer variable size" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  // Each rank sends (rank+1) elements
  std::vector<int> local(rank + 1);
  for (int i = 0; i < rank + 1; ++i) {
    local[i] = rank * 100 + i;
  }
  auto collected = mpi.collectContainer(local);

  if (rank == 0) {
    REQUIRE(collected.size() == static_cast<std::size_t>(size));
    for (int r = 0; r < size; ++r) {
      REQUIRE(collected[r].size() == static_cast<std::size_t>(r + 1));
      for (int i = 0; i < r + 1; ++i) {
        CHECK(collected[r][i] == r * 100 + i);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// broadcast<T>
// ---------------------------------------------------------------------------

TEST_CASE("MPI broadcast double" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();

  double value = (rank == 0) ? 3.14159 : 0.0;
  mpi.broadcast(&value, 0);

  CHECK(value == doctest::Approx(3.14159));
}

TEST_CASE("MPI broadcast int" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();

  int value = (rank == 0) ? 42 : -1;
  mpi.broadcast(&value, 0);

  CHECK(value == 42);
}

TEST_CASE("MPI broadcast from non-root" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  if (size < 2) {
    return; // Need at least 2 ranks
  }

  // Broadcast from rank 1
  int value = (rank == 1) ? 999 : 0;
  mpi.broadcast(&value, 1);

  CHECK(value == 999);
}

// ---------------------------------------------------------------------------
// broadcastContainer
// ---------------------------------------------------------------------------

TEST_CASE("MPI broadcastContainer" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();

  std::vector<double> data;
  if (rank == 0) {
    data = {1.1, 2.2, 3.3, 4.4, 5.5};
  }

  mpi.broadcastContainer(data, 0);

  // All ranks should now have the data
  REQUIRE(data.size() == 5);
  CHECK(data[0] == doctest::Approx(1.1));
  CHECK(data[1] == doctest::Approx(2.2));
  CHECK(data[2] == doctest::Approx(3.3));
  CHECK(data[3] == doctest::Approx(4.4));
  CHECK(data[4] == doctest::Approx(5.5));
}

TEST_CASE("MPI broadcastContainer empty from root" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();

  std::vector<int> data;
  if (rank == 0) {
    data = {};
  } else {
    data = {999}; // Non-root has garbage
  }

  mpi.broadcastContainer(data, 0);

  CHECK(data.empty());
}

// ---------------------------------------------------------------------------
// MPI rank/size consistency
// ---------------------------------------------------------------------------

TEST_CASE("MPI rank and size consistency" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;

  CHECK(mpi.rank() >= 0);
  CHECK(mpi.rank() < mpi.size());
  CHECK(mpi.size() >= 1);
}

// ---------------------------------------------------------------------------
// barrier does not deadlock
// ---------------------------------------------------------------------------

TEST_CASE("MPI barrier" * doctest::test_suite("mpi")) {
  seissol::Mpi::barrier(MPI_COMM_WORLD);
  CHECK(true);
}

} // namespace seissol::unit_test
