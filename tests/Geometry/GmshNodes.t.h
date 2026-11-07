// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Geometry/GmshNodes.h"

#include <cstddef>
#include <vector>

namespace seissol::unit_test::gmshnodes {

using seissol::geometry::LatticePoint;

/// the nodes of the tetrahedra of orders one to five as gmsh 4.15.2 numbers them
/// (gmsh.model.mesh.getElementProperties), times the order
inline auto gmshReference() -> std::vector<std::vector<LatticePoint>> {
  return {
      // order 1
      {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
      // order 2
      {{0, 0, 0},
       {2, 0, 0},
       {0, 2, 0},
       {0, 0, 2},
       {1, 0, 0},
       {1, 1, 0},
       {0, 1, 0},
       {0, 0, 1},
       {0, 1, 1},
       {1, 0, 1}},
      // order 3
      {{0, 0, 0}, {3, 0, 0}, {0, 3, 0}, {0, 0, 3}, {1, 0, 0}, {2, 0, 0}, {2, 1, 0},
       {1, 2, 0}, {0, 2, 0}, {0, 1, 0}, {0, 0, 2}, {0, 0, 1}, {0, 1, 2}, {0, 2, 1},
       {1, 0, 2}, {2, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}},
      // order 4
      {{0, 0, 0}, {4, 0, 0}, {0, 4, 0}, {0, 0, 4}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0},
       {3, 1, 0}, {2, 2, 0}, {1, 3, 0}, {0, 3, 0}, {0, 2, 0}, {0, 1, 0}, {0, 0, 3},
       {0, 0, 2}, {0, 0, 1}, {0, 1, 3}, {0, 2, 2}, {0, 3, 1}, {1, 0, 3}, {2, 0, 2},
       {3, 0, 1}, {1, 1, 0}, {1, 2, 0}, {2, 1, 0}, {1, 0, 1}, {2, 0, 1}, {1, 0, 2},
       {0, 1, 1}, {0, 1, 2}, {0, 2, 1}, {1, 1, 2}, {2, 1, 1}, {1, 2, 1}, {1, 1, 1}},
      // order 5
      {{0, 0, 0}, {5, 0, 0}, {0, 5, 0}, {0, 0, 5}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {4, 0, 0},
       {4, 1, 0}, {3, 2, 0}, {2, 3, 0}, {1, 4, 0}, {0, 4, 0}, {0, 3, 0}, {0, 2, 0}, {0, 1, 0},
       {0, 0, 4}, {0, 0, 3}, {0, 0, 2}, {0, 0, 1}, {0, 1, 4}, {0, 2, 3}, {0, 3, 2}, {0, 4, 1},
       {1, 0, 4}, {2, 0, 3}, {3, 0, 2}, {4, 0, 1}, {1, 1, 0}, {1, 3, 0}, {3, 1, 0}, {1, 2, 0},
       {2, 2, 0}, {2, 1, 0}, {1, 0, 1}, {3, 0, 1}, {1, 0, 3}, {2, 0, 1}, {2, 0, 2}, {1, 0, 2},
       {0, 1, 1}, {0, 1, 3}, {0, 3, 1}, {0, 1, 2}, {0, 2, 2}, {0, 2, 1}, {1, 1, 3}, {3, 1, 1},
       {1, 3, 1}, {2, 1, 2}, {2, 2, 1}, {1, 2, 2}, {1, 1, 1}, {2, 1, 1}, {1, 2, 1}, {1, 1, 2}},
  };
}

TEST_CASE("The nodes of a tetrahedron are numbered as gmsh numbers them") {
  const auto reference = gmshReference();
  for (std::size_t order = 1; order <= reference.size(); ++order) {
    CAPTURE(order);
    const auto lattice = seissol::geometry::gmshTetrahedronLattice(order);
    REQUIRE(lattice.size() == seissol::geometry::lagrangeNodeCount(order));
    REQUIRE(lattice == reference[order - 1]);
  }
}

} // namespace seissol::unit_test::gmshnodes
