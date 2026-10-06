// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_GMSHNODES_H_
#define SEISSOL_SRC_GEOMETRY_GMSHNODES_H_

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::geometry {

/// A node of the lattice of a tetrahedron of order p: its reference coordinates times p.
using LatticePoint = std::array<std::size_t, Cell::Dim>;

/**
 * The nodes of a Lagrange tetrahedron of `order` in the order gmsh numbers them, which is the order
 * a mesh file of PUMGen gives them in: the four vertices, then the nodes inside the edges, edge by
 * edge from the first vertex of the edge to its second, then the nodes inside the faces, face by
 * face, and last the ones inside the cell. The nodes inside a face are numbered as gmsh numbers a
 * triangle of order - 3 put into the face, and the ones inside the cell as gmsh numbers a
 * tetrahedron of order - 4 put into it.
 */
auto gmshTetrahedronLattice(std::size_t order) -> std::vector<LatticePoint>;

/// The number of nodes of a Lagrange tetrahedron of `order`.
constexpr auto lagrangeNodeCount(std::size_t order) -> std::size_t {
  return (order + 1) * (order + 2) * (order + 3) / 6;
}

/**
 * Where each node of IsoparametricTransform::latticeNodes(order) stands in the node list of gmsh,
 * for a cell whose local vertex k is the vertex vertexOrder[k] of the cell as the file gives it
 * (the convention of geometry::VertexOrder). The nodes of the file are positioned relative to the
 * vertices of the file, so a cell whose vertices the reader puts into another order needs its
 * other nodes put into the corresponding order as well.
 */
auto gmshToLattice(std::size_t order,
                   const std::array<std::uint8_t, Cell::NumVertices>& vertexOrder)
    -> std::vector<std::size_t>;

/**
 * The nodes of a cell on the lattice of `targetOrder`, in the order of
 * IsoparametricTransform::latticeNodes(targetOrder), from what a mesh file gives for it: its
 * vertices in the vertex order of the file, and its other `order` nodes in the order of gmsh, as
 * three coordinates each. A cell of an order below the target order is the same map written with
 * more nodes, which is how the cells of a mesh whose cells differ in order get one order. Its local
 * vertex k is the vertex vertexOrder[k] of the file.
 */
auto latticeNodesFromGmsh(
    std::size_t order,
    std::size_t targetOrder,
    const std::array<CellTransform::VectorEigenT, Cell::NumVertices>& vertices,
    const double* otherNodes,
    const std::array<std::uint8_t, Cell::NumVertices>& vertexOrder)
    -> std::vector<CellTransform::VectorEigenT>;

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_GMSHNODES_H_
