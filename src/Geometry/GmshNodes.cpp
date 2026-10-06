// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GmshNodes.h"

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/IsoparametricTransform.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <utils/logger.h>
#include <vector>

namespace seissol::geometry {

namespace {

/// A node by its barycentric coordinates times the order, for a simplex of N vertices
template <std::size_t N>
using Barycentric = std::array<std::size_t, N>;

/**
 * The nodes of a Lagrange simplex of N vertices and `order` in the order of gmsh: the vertices,
 * the nodes inside the edges in the order of `edges`, each from its first vertex to its second,
 * the nodes inside the faces of a tetrahedron, and those inside the simplex itself, which are the
 * nodes of the simplex of order - N lifted by one in every coordinate.
 */
template <std::size_t N>
auto gmshSimplex(std::size_t order,
                 const std::vector<std::array<std::size_t, 2>>& edges,
                 const std::vector<std::array<std::size_t, 3>>& faces)
    -> std::vector<Barycentric<N>> {
  std::vector<Barycentric<N>> nodes;
  if (order == 0) {
    nodes.push_back(Barycentric<N>{});
    return nodes;
  }
  for (std::size_t vertex = 0; vertex < N; ++vertex) {
    Barycentric<N> node{};
    node[vertex] = order;
    nodes.push_back(node);
  }
  for (const auto& edge : edges) {
    for (std::size_t i = 1; i < order; ++i) {
      Barycentric<N> node{};
      node[edge[0]] = order - i;
      node[edge[1]] = i;
      nodes.push_back(node);
    }
  }
  if (order >= 3) {
    // the nodes inside a face are the ones of a triangle of order - 3, lifted by one, with the
    // vertices of the triangle put onto the vertices of the face in the order the face lists them
    for (const auto& face : faces) {
      for (const auto& inner : gmshSimplex<3>(order - 3, {{0, 1}, {1, 2}, {2, 0}}, {})) {
        Barycentric<N> node{};
        for (std::size_t k = 0; k < 3; ++k) {
          node[face[k]] = inner[k] + 1;
        }
        nodes.push_back(node);
      }
    }
  }
  if (order >= N) {
    for (const auto& inner : gmshSimplex<N>(order - N, edges, faces)) {
      Barycentric<N> node{};
      for (std::size_t k = 0; k < N; ++k) {
        node[k] = inner[k] + 1;
      }
      nodes.push_back(node);
    }
  }
  return nodes;
}

// the edges and faces of a tetrahedron as gmsh numbers them
const std::vector<std::array<std::size_t, 2>> GmshEdges = {
    {0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}};
const std::vector<std::array<std::size_t, 3>> GmshFaces = {
    {0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}};

} // namespace

auto gmshTetrahedronLattice(std::size_t order) -> std::vector<LatticePoint> {
  std::vector<LatticePoint> lattice;
  for (const auto& node : gmshSimplex<Cell::NumVertices>(order, GmshEdges, GmshFaces)) {
    lattice.push_back({node[1], node[2], node[3]});
  }
  return lattice;
}

auto gmshToLattice(std::size_t order,
                   const std::array<std::uint8_t, Cell::NumVertices>& vertexOrder)
    -> std::vector<std::size_t> {
  std::map<LatticePoint, std::size_t> gmshIndex;
  const auto gmsh = gmshTetrahedronLattice(order);
  for (std::size_t i = 0; i < gmsh.size(); ++i) {
    gmshIndex[gmsh[i]] = i;
  }

  const auto lattice = IsoparametricTransform::latticeNodes(order);
  std::vector<std::size_t> indices(lattice.size());
  for (std::size_t i = 0; i < lattice.size(); ++i) {
    // the node by its barycentric coordinates on the local vertices, then on those of the file
    std::array<std::size_t, Cell::NumVertices> local{};
    std::size_t sum = 0;
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      local[d + 1] = static_cast<std::size_t>(std::lround(lattice[i](d) * order));
      sum += local[d + 1];
    }
    local[0] = order - sum;
    std::array<std::size_t, Cell::NumVertices> file{};
    for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
      file[vertexOrder[k]] = local[k];
    }
    const auto found = gmshIndex.find({file[1], file[2], file[3]});
    if (found == gmshIndex.end()) {
      logError() << "A node of the lattice of order" << order << "has no node of gmsh.";
    }
    indices[i] = found->second;
  }
  return indices;
}

auto latticeNodesFromGmsh(
    std::size_t order,
    std::size_t targetOrder,
    const std::array<CellTransform::VectorEigenT, Cell::NumVertices>& vertices,
    const double* otherNodes,
    const std::array<std::uint8_t, Cell::NumVertices>& vertexOrder)
    -> std::vector<CellTransform::VectorEigenT> {
  if (order < 1 || order > targetOrder) {
    logError() << "A cell of order" << order << "cannot be written on the lattice of order"
               << targetOrder << ".";
  }
  const auto indices = gmshToLattice(order, vertexOrder);
  std::vector<CellTransform::VectorEigenT> nodes(indices.size());
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const auto index = indices[i];
    if (index < Cell::NumVertices) {
      nodes[i] = vertices[index];
    } else {
      const auto* node = otherNodes + 3 * (index - Cell::NumVertices);
      nodes[i] = CellTransform::VectorEigenT(node[0], node[1], node[2]);
    }
  }
  if (order == targetOrder) {
    return nodes;
  }
  // the map of the cell, written with the nodes of the target order
  const IsoparametricTransform transform(order, nodes);
  const auto lattice = IsoparametricTransform::latticeNodes(targetOrder);
  std::vector<CellTransform::VectorEigenT> elevated(lattice.size());
  for (std::size_t i = 0; i < lattice.size(); ++i) {
    elevated[i] = transform.refToSpace(lattice[i]);
  }
  return elevated;
}

} // namespace seissol::geometry
