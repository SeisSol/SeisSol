// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GmshNodes.h"

#include "Common/Constants.h"

#include <array>
#include <cstddef>
#include <vector>

namespace seissol::geometry {

namespace {

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

} // namespace seissol::geometry
