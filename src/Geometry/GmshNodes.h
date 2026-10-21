// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_GEOMETRY_GMSHNODES_H_
#define SEISSOL_SRC_GEOMETRY_GMSHNODES_H_

#include "Common/Constants.h"

#include <array>
#include <cstddef>
#include <vector>

namespace seissol::geometry {

using LatticePoint = std::array<std::size_t, Cell::Dim>;

/**
 * The nodes of a Lagrange tetrahedron of `order` in the order gmsh numbers them, which is the order
 * a mesh file of PUMGen gives them in: the four vertices, then the nodes inside the edges, edge by
 * edge from the first vertex of the edge to its second, then the nodes inside the faces, face by
 * face, and last the ones inside the cell. The nodes inside a face are numbered as gmsh numbers a
 * triangle of order - 3 put into the face, and the ones inside the cell as gmsh numbers a
 * tetrahedron of order - 4 put into it.
 *
 * A node is given by its coordinates on the reference cell times the order.
 */
auto gmshTetrahedronLattice(std::size_t order) -> std::vector<LatticePoint>;

/// the nodes of a Lagrange tetrahedron of `order`, vertices included
constexpr auto lagrangeNodeCount(std::size_t order) -> std::size_t {
  return (order + 1) * (order + 2) * (order + 3) / 6;
}

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_GMSHNODES_H_
