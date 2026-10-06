// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "IsoparametricTransform.h"

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::geometry {

namespace {

using VectorEigenT = CellTransform::VectorEigenT;

auto power(double base, int exponent) -> double {
  double result = 1.0;
  for (int i = 0; i < exponent; ++i) {
    result *= base;
  }
  return result;
}

auto monomial(const std::array<int, Cell::Dim>& exponent, const VectorEigenT& point) -> double {
  double result = 1.0;
  for (std::size_t d = 0; d < Cell::Dim; ++d) {
    result *= power(point(d), exponent[d]);
  }
  return result;
}

auto monomialDerivative(const std::array<int, Cell::Dim>& exponent,
                        const VectorEigenT& point,
                        std::size_t direction) -> double {
  if (exponent[direction] == 0) {
    return 0.0;
  }
  double result = exponent[direction];
  for (std::size_t d = 0; d < Cell::Dim; ++d) {
    result *= power(point(d), exponent[d] - (d == direction ? 1 : 0));
  }
  return result;
}

/// the barycentric coordinates of a reference point, the first one belonging to vertex zero
auto barycentric(const VectorEigenT& point) -> std::array<double, Cell::NumVertices> {
  return {1.0 - point(0) - point(1) - point(2), point(0), point(1), point(2)};
}

} // namespace

auto IsoparametricTransform::latticeNodes(std::size_t order) -> std::vector<VectorEigenT> {
  if (order == 0) {
    logError() << "An isoparametric cell needs at least the order one.";
  }
  std::vector<VectorEigenT> nodes{
      VectorEigenT(0, 0, 0), VectorEigenT(1, 0, 0), VectorEigenT(0, 1, 0), VectorEigenT(0, 0, 1)};
  const auto scale = static_cast<double>(order);
  const auto top = static_cast<int>(order);
  for (int k = 0; k <= top; ++k) {
    for (int j = 0; j + k <= top; ++j) {
      for (int i = 0; i + j + k <= top; ++i) {
        const bool vertex = (i == 0 && j == 0 && k == 0) || (i == top) || (j == top) || (k == top);
        if (!vertex) {
          nodes.emplace_back(i / scale, j / scale, k / scale);
        }
      }
    }
  }
  return nodes;
}

IsoparametricTransform::IsoparametricTransform(std::size_t order,
                                               const std::vector<VectorEigenT>& nodes)
    : order_(order), nodes_(nodes) {
  const auto reference = latticeNodes(order);
  if (nodes.size() != reference.size()) {
    logError() << "An isoparametric cell of order" << order << "has" << reference.size()
               << "nodes, but" << nodes.size() << "were given.";
  }
  const auto top = static_cast<int>(order);
  for (int k = 0; k <= top; ++k) {
    for (int j = 0; j + k <= top; ++j) {
      for (int i = 0; i + j + k <= top; ++i) {
        exponents_.push_back({i, j, k});
      }
    }
  }

  const auto count = static_cast<Eigen::Index>(reference.size());
  Eigen::MatrixXd vandermonde(count, count);
  Eigen::MatrixXd shifted(count, static_cast<Eigen::Index>(Cell::Dim));
  origin_ = nodes[0];
  for (Eigen::Index row = 0; row < count; ++row) {
    for (Eigen::Index column = 0; column < count; ++column) {
      vandermonde(row, column) = monomial(exponents_[column], reference[row]);
    }
    shifted.row(row) = (nodes[row] - origin_).transpose();
  }
  coefficients_ = vandermonde.fullPivLu().solve(shifted).transpose();
}

auto IsoparametricTransform::fromVertices(
    const std::array<VectorEigenT, Cell::NumVertices>& vertices, std::size_t order)
    -> IsoparametricTransform {
  const auto reference = latticeNodes(order);
  std::vector<VectorEigenT> nodes;
  nodes.reserve(reference.size());
  for (const auto& point : reference) {
    const auto weights = barycentric(point);
    VectorEigenT node = VectorEigenT::Zero();
    for (std::size_t vertex = 0; vertex < Cell::NumVertices; ++vertex) {
      node += weights[vertex] * vertices[vertex];
    }
    nodes.push_back(node);
  }
  return {order, nodes};
}

auto IsoparametricTransform::fromEdgeMidpoints(
    const std::array<VectorEigenT, Cell::NumVertices>& vertices,
    const std::array<VectorEigenT, 6>& midpoints) -> IsoparametricTransform {
  // the edges in the order the midpoints are given in
  constexpr std::array<std::array<std::size_t, 2>, 6> Edges{
      {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
  const auto reference = latticeNodes(2);
  std::vector<VectorEigenT> nodes;
  nodes.reserve(reference.size());
  for (const auto& point : reference) {
    const auto weights = barycentric(point);
    std::array<std::size_t, 2> touched{};
    std::size_t count = 0;
    for (std::size_t vertex = 0; vertex < Cell::NumVertices; ++vertex) {
      if (std::abs(weights[vertex]) > 1e-12 && count < 2) {
        touched[count] = vertex;
      }
      count += std::abs(weights[vertex]) > 1e-12 ? 1 : 0;
    }
    if (count == 1) {
      nodes.push_back(vertices[touched[0]]);
      continue;
    }
    for (std::size_t edge = 0; edge < Edges.size(); ++edge) {
      if (Edges[edge] == touched) {
        nodes.push_back(midpoints[edge]);
      }
    }
  }
  return {2, nodes};
}

auto IsoparametricTransform::refToSpace(const VectorEigenT& input) const -> VectorEigenT {
  VectorEigenT result = origin_;
  for (std::size_t i = 0; i < exponents_.size(); ++i) {
    result += monomial(exponents_[i], input) * coefficients_.col(static_cast<Eigen::Index>(i));
  }
  return result;
}

auto IsoparametricTransform::refToSpaceJacobian(const VectorEigenT& input) const -> MatrixEigenT {
  MatrixEigenT jacobian = MatrixEigenT::Zero();
  for (std::size_t i = 0; i < exponents_.size(); ++i) {
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      jacobian.col(static_cast<Eigen::Index>(d)) += monomialDerivative(exponents_[i], input, d) *
                                                    coefficients_.col(static_cast<Eigen::Index>(i));
    }
  }
  return jacobian;
}

IsoparametricFaceTransform::IsoparametricFaceTransform(const IsoparametricTransform& cell,
                                                       const ReferenceFaceMap& embedding)
    : map_(embedding), cell_(cell) {}

auto IsoparametricFaceTransform::refToSpace(const FaceVectorT& input) const -> VectorT {
  return cell_.refToSpace(map_.faceToCell(input));
}

auto IsoparametricFaceTransform::refToCell(const FaceVectorT& input) const -> VectorT {
  return map_.faceToCell(input);
}

auto IsoparametricFaceTransform::refToSpaceJacobian(const FaceVectorT& input) const -> JacobianT {
  return cell_.refToSpaceJacobian(map_.faceToCell(input)) * map_.embedding();
}

} // namespace seissol::geometry
