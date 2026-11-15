// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/MPI.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <string>
#include <utility>
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
// allreduce, allreduceContainer
// ---------------------------------------------------------------------------

TEST_CASE("MPI allreduce of a value" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  CHECK(mpi.allreduce(rank, MPI_SUM) == size * (size - 1) / 2);
  CHECK(mpi.allreduce(static_cast<std::size_t>(rank), MPI_MAX) ==
        static_cast<std::size_t>(size - 1));
  CHECK(mpi.allreduce(rank + 0.5, MPI_MIN) == doctest::Approx(0.5));
  CHECK(mpi.allreduce(rank == 0, MPI_LOR));
  CHECK_FALSE(mpi.allreduce(rank == 0 && size > 1, MPI_LAND));
}

TEST_CASE("MPI allreduceContainer reduces in place" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  std::vector<double> values{1.0, static_cast<double>(rank)};
  mpi.allreduceContainer(values, MPI_SUM);
  CHECK(values[0] == doctest::Approx(size));
  CHECK(values[1] == doctest::Approx(size * (size - 1) / 2.0));

  std::array<std::uint64_t, 2> extremes{static_cast<std::uint64_t>(rank),
                                        static_cast<std::uint64_t>(rank)};
  mpi.allreduceContainer(extremes, MPI_MAX);
  CHECK(extremes[0] == static_cast<std::uint64_t>(size - 1));
  CHECK(extremes[1] == static_cast<std::uint64_t>(size - 1));
}

TEST_CASE("MPI allreduceContainer finds where the extremes are" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  static_assert(seissol::Mpi::pairMatchesMpiLayout<double, int>());
  // the smallest value lies on the last rank, the largest on the first
  std::vector<std::pair<double, int>> smallest{{static_cast<double>(size - rank), rank}};
  std::vector<std::pair<double, int>> largest{{static_cast<double>(size - rank), rank}};
  mpi.allreduceContainer(smallest, MPI_MINLOC);
  mpi.allreduceContainer(largest, MPI_MAXLOC);
  CHECK(smallest[0].first == doctest::Approx(1.0));
  CHECK(smallest[0].second == size - 1);
  CHECK(largest[0].first == doctest::Approx(size));
  CHECK(largest[0].second == 0);
}

// ---------------------------------------------------------------------------
// scan
// ---------------------------------------------------------------------------

TEST_CASE("MPI scan" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();

  const auto value = static_cast<std::size_t>(rank + 1);
  const auto rankSize = static_cast<std::size_t>(rank);
  CHECK(mpi.scan(value, MPI_SUM, true) == (rankSize + 1) * (rankSize + 2) / 2);
  // defined on the first rank as well, where MPI_Exscan leaves it undefined
  CHECK(mpi.scan(value, MPI_SUM, false) == rankSize * (rankSize + 1) / 2);
}

TEST_CASE("MPI scanContainer" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const auto rank = static_cast<std::size_t>(mpi.rank());

  std::vector<std::size_t> inclusive{1, rank};
  auto exclusive = inclusive;
  mpi.scanContainer(inclusive, MPI_SUM, true);
  mpi.scanContainer(exclusive, MPI_SUM, false);
  CHECK(inclusive[0] == rank + 1);
  CHECK(inclusive[1] == rank * (rank + 1) / 2);
  CHECK(exclusive[0] == rank);
  CHECK(exclusive[1] == rank * (rank - 1) / 2);
}

// ---------------------------------------------------------------------------
// reduce, reduceContainer, allgather
// ---------------------------------------------------------------------------

TEST_CASE("MPI reduce to one rank" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();
  const int root = size - 1;

  const auto sum = mpi.reduce(rank, MPI_SUM, root);
  std::vector<double> values{1.0, static_cast<double>(rank)};
  mpi.reduceContainer(values, MPI_MAX, root);
  if (rank == root) {
    CHECK(sum == size * (size - 1) / 2);
    CHECK(values[1] == doctest::Approx(size - 1));
  } else {
    // the others keep their own values
    CHECK(sum == 0);
    CHECK(values[1] == doctest::Approx(rank));
  }
  CHECK(values[0] == doctest::Approx(1.0));
}

TEST_CASE("MPI allgather" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;

  const auto values = mpi.allgather(10 * mpi.rank());
  REQUIRE(values.size() == static_cast<std::size_t>(mpi.size()));
  for (int other = 0; other < mpi.size(); ++other) {
    CHECK(values[other] == 10 * other);
  }
}

// ---------------------------------------------------------------------------
// other communicators
// ---------------------------------------------------------------------------

namespace {
/// the ranks of MPI_COMM_WORLD split by their parity, the even ones first
MPI_Comm splitByParity() {
  MPI_Comm comm = MPI_COMM_NULL;
  const int rank = seissol::Mpi::rankOf(MPI_COMM_WORLD);
  MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &comm);
  return comm;
}
} // namespace

TEST_CASE("MPI collectives on a communicator of their own" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int worldRank = mpi.rank();
  const int worldSize = mpi.size();
  MPI_Comm comm = splitByParity();

  const int rank = seissol::Mpi::rankOf(comm);
  const int size = seissol::Mpi::sizeOf(comm);
  CHECK(rank == worldRank / 2);
  CHECK(size == (worldSize + 1 - worldRank % 2) / 2);

  // the gathers have as many entries as the communicator has ranks, on its own first rank
  const auto ranks = mpi.collect(worldRank, comm);
  const auto names = mpi.collectContainer(std::vector<int>(rank + 1, worldRank), comm);
  if (rank == 0) {
    REQUIRE(ranks.size() == static_cast<std::size_t>(size));
    REQUIRE(names.size() == static_cast<std::size_t>(size));
    for (int other = 0; other < size; ++other) {
      CHECK(ranks[other] == 2 * other + worldRank % 2);
      CHECK(names[other] == std::vector<int>(other + 1, 2 * other + worldRank % 2));
    }
  }

  // from the last rank of the communicator, which is not the first one of the world
  std::vector<int> data;
  if (rank == size - 1) {
    data = {worldRank, worldRank};
  }
  mpi.broadcastContainer(data, size - 1, comm);
  const int lastWorldRank = 2 * (size - 1) + worldRank % 2;
  CHECK(data == std::vector<int>{lastWorldRank, lastWorldRank});

  CHECK(mpi.allreduce(1, MPI_SUM, comm) == size);
  CHECK(mpi.scan(1, MPI_SUM, false, comm) == rank);
  mpi.barrier(comm);

  MPI_Comm_free(&comm);
}

TEST_CASE("MPI collectives take the communicator of the thread" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  MPI_Comm comm = splitByParity();

  seissol::Mpi::setThreadComm(comm);
  CHECK(mpi.comm() == comm);
  const int size = seissol::Mpi::sizeOf(comm);
  CHECK(mpi.allreduce(1, MPI_SUM) == size);
  const auto ranks = mpi.collect(mpi.rank());
  if (seissol::Mpi::rankOf(comm) == 0) {
    CHECK(ranks.size() == static_cast<std::size_t>(size));
  }
  seissol::Mpi::setThreadComm(MPI_COMM_NULL);
  CHECK(mpi.allreduce(1, MPI_SUM) == mpi.size());

  MPI_Comm_free(&comm);
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
  seissol::Mpi::mpi.barrier(MPI_COMM_WORLD);
  CHECK(true);
}

} // namespace seissol::unit_test
