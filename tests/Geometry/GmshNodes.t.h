// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/GmshNodes.h"
#include "Geometry/IsoparametricTransform.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

namespace seissol::unit_test::gmshnodes {

using seissol::geometry::CellTransform;
using seissol::geometry::IsoparametricTransform;
using seissol::geometry::LatticePoint;
using VectorT = CellTransform::VectorEigenT;

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

/// a map of the reference cell of degree `order` in every component
inline auto polynomialMap(std::size_t order, std::mt19937& rng) {
  std::normal_distribution<double> gauss(0.0, 1.0);
  Eigen::Matrix3d linear;
  linear << 2.0, 0.3, -0.2, 0.1, 1.7, 0.4, -0.3, 0.2, 2.2;
  std::vector<std::array<double, 4>> terms;
  for (std::size_t i = 0; i < 3 * order; ++i) {
    terms.push_back({0.05 * gauss(rng), 0.05 * gauss(rng), 0.05 * gauss(rng), 0.0});
  }
  return [=](const VectorT& xi) -> VectorT {
    VectorT x = VectorT(0.3, -1.2, 5.0) + linear * xi;
    if (order >= 2) {
      for (std::size_t i = 0; i < terms.size(); ++i) {
        // monomials of degree two up to the order, one per term
        const std::size_t degree = 2 + i % (order - 1);
        const double value = std::pow(xi((i + 0) % 3), degree - 1) * xi((i + 1) % 3);
        x += value * VectorT(terms[i][0], terms[i][1], terms[i][2]);
      }
    }
    return x;
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

TEST_CASE("The nodes of gmsh give the map of the cell in any vertex order") {
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937 rng(20261006);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (std::size_t order = 1; order <= 4; ++order) {
    for (std::size_t target = order; target <= 4; ++target) {
      CAPTURE(order);
      CAPTURE(target);
      const auto map = polynomialMap(order, rng);
      // the cell as the file gives it: its vertices and its other nodes, in the order of gmsh
      const auto gmsh = seissol::geometry::gmshTetrahedronLattice(order);
      std::array<VectorT, Cell::NumVertices> vertices;
      std::vector<double> otherNodes;
      for (std::size_t i = 0; i < gmsh.size(); ++i) {
        const VectorT xi = VectorT(gmsh[i][0], gmsh[i][1], gmsh[i][2]) / static_cast<double>(order);
        const VectorT x = map(xi);
        if (i < Cell::NumVertices) {
          vertices[i] = x;
        } else {
          otherNodes.insert(otherNodes.end(), {x(0), x(1), x(2)});
        }
      }

      std::array<std::size_t, Cell::NumVertices> permutation{0, 1, 2, 3};
      do {
        std::array<std::uint8_t, Cell::NumVertices> vertexOrder{};
        for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
          vertexOrder[k] = static_cast<std::uint8_t>(permutation[k]);
        }
        const auto nodes = seissol::geometry::latticeNodesFromGmsh(
            order, target, vertices, otherNodes.data(), vertexOrder);
        REQUIRE(nodes.size() == seissol::geometry::lagrangeNodeCount(target));
        for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
          REQUIRE((nodes[k] - vertices[vertexOrder[k]]).norm() == AbsApprox(0.0).epsilon(1e-12));
        }
        const IsoparametricTransform cell(target, nodes);
        for (int sample = 0; sample < 8; ++sample) {
          // a point of the cell by its barycentric coordinates on the local vertices, and on the
          // vertices of the file
          std::array<double, Cell::NumVertices> local{};
          double sum = 0;
          for (auto& weight : local) {
            weight = -std::log(unit(rng) + 1e-300);
            sum += weight;
          }
          std::array<double, Cell::NumVertices> file{};
          for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
            local[k] /= sum;
            file[vertexOrder[k]] = local[k];
          }
          const VectorT expected = map(VectorT(file[1], file[2], file[3]));
          const VectorT actual = cell.refToSpace(VectorT(local[1], local[2], local[3]));
          REQUIRE((actual - expected).norm() == AbsApprox(0.0).epsilon(1e-10));
        }
      } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
  }
}

} // namespace seissol::unit_test::gmshnodes
