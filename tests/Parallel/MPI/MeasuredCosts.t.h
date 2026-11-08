// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/Clustering/VertexWeights/MeasuredCellCosts.h"
#include "Initializer/Clustering/VertexWeights/MeasuredWeights.h"
#include "Parallel/MPI.h"
#include "TestHelper.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace measuredcoststest {

constexpr std::size_t CellsPerRank = 7;

/**
 * The costs a run measured, as a partition leaves them: rank r holds the cells r, r + size, ...
 * (backwards, so that no rank holds its cells in order), each costing a tenth of its id.
 */
inline initializer::MeasuredCellCosts scattered() {
  const auto rank = static_cast<std::size_t>(Mpi::mpi.rank());
  const auto size = static_cast<std::size_t>(Mpi::mpi.size());
  initializer::MeasuredCellCosts costs;
  for (std::size_t cell = rank; cell < CellsPerRank * size; cell += size) {
    costs.ids.push_back(cell);
    costs.costs.push_back(0.1 * static_cast<double>(cell));
  }
  std::ranges::reverse(costs.ids);
  std::ranges::reverse(costs.costs);
  return costs;
}

/// the run of cells the rank reads from the mesh file: ranks before it read one cell fewer each
inline std::pair<std::size_t, std::size_t> runOfRank() {
  const auto rank = static_cast<std::size_t>(Mpi::mpi.rank());
  const auto size = static_cast<std::size_t>(Mpi::mpi.size());
  const auto total = CellsPerRank * size;
  // CellsPerRank - 1 for all but the last rank, which reads the rest
  const auto first = rank * (CellsPerRank - 1);
  const auto count = rank + 1 == size ? total - first : CellsPerRank - 1;
  return {first, count};
}

} // namespace measuredcoststest

TEST_CASE("Measured costs reach the ranks that read the cells" * doctest::test_suite("mpi")) {
  using namespace measuredcoststest;
  const auto costs = scattered();
  const auto [first, count] = runOfRank();

  const auto received = initializer::distributeCellCosts(costs, first, count, Mpi::mpi.comm());
  REQUIRE(received.size() == count);
  for (std::size_t i = 0; i < count; ++i) {
    CHECK(received[i] == doctest::Approx(0.1 * static_cast<double>(first + i)));
  }
}

TEST_CASE("Measured costs go through a file" * doctest::test_suite("mpi")) {
  using namespace measuredcoststest;

  std::string dir;
  if (Mpi::mpi.rank() == 0) {
    dir = createTemporaryDirectory("seissolcosts");
  }
  Mpi::mpi.broadcastContainer(dir, 0);
  const auto file = dir + "/costs.h5";

  const auto costs = scattered();
  initializer::writeCellCosts(file, costs, 2.5, Mpi::mpi.comm());

  const auto [first, count] = runOfRank();
  const auto read = initializer::readCellCosts(file, first, count, Mpi::mpi.comm());
  REQUIRE(read.size() == count);
  for (std::size_t i = 0; i < count; ++i) {
    CHECK(read[i] == doctest::Approx(0.1 * static_cast<double>(first + i)));
  }

  MPI_Barrier(Mpi::mpi.comm());
  if (Mpi::mpi.rank() == 0) {
    std::filesystem::remove_all(dir);
  }
}

TEST_CASE("The time clusters go to groups with a share of the costs each" *
          doctest::test_suite("mpi")) {
  using initializer::constraintGroups;
  using Groups = std::vector<std::size_t>;
  // two tiny clusters with small time steps join the next one, the bulk is a group of its own
  CHECK(constraintGroups({1e-5, 0.002, 0.2, 0.798}, 0.1) == Groups{0, 0, 0, 1});
  // a last group with too little of the costs joins the one before
  CHECK(constraintGroups({0.5, 0.45, 0.05}, 0.1) == Groups{0, 1, 1});
  CHECK(constraintGroups({0.25, 0.25, 0.25, 0.25}, 0.1) == Groups{0, 1, 2, 3});
  // a share above all costs: one group
  CHECK(constraintGroups({0.3, 0.7}, 2.0) == Groups{0, 0});
  CHECK(constraintGroups({1.0}, 0.1) == Groups{0});
}

} // namespace seissol::unit_test
