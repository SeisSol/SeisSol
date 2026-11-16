// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/InitProcedure/Internal/NeighborLinks.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Parallel/MPI.h"

#include <cmath>
#include <set>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("The links to the neighbor ranks of a halo" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  // a ring: each rank exchanges cells with the ranks before and after it
  std::set<int> neighbors;
  initializer::MeshLayout layout(1);
  if (size > 1) {
    for (const int neighbor : {(rank + size - 1) % size, (rank + 1) % size}) {
      layout[0].regions.emplace_back(0, 0, 0, 1, neighbor);
      neighbors.insert(neighbor);
    }
  }

  const auto links = initializer::internal::measureNeighborLinks(layout);
  REQUIRE(links.size() == neighbors.size());
  for (const auto neighbor : neighbors) {
    REQUIRE(links.contains(neighbor));
    const auto& link = links.at(neighbor);
    CHECK(link.latency > 0);
    CHECK(std::isfinite(link.latency));
    CHECK(link.bytesPerSecond > 0);
    // a larger message takes longer, or as long
    CHECK(link.time(1 << 20) >= link.time(8));
  }
}

} // namespace seissol::unit_test
