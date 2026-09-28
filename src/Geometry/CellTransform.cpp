// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CellTransform.h"

#include "Common/Constants.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"

#include <Eigen/Dense>
#include <Eigen/LU>
#include <array>
#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::geometry {

CellTransform::~CellTransform() = default;

auto CellTransform::refToSpace(const VectorT& input) const -> VectorT {
  const auto outputEigen = refToSpace(VectorEigenT(input.data()));
  VectorT output{};
  for (std::size_t i = 0; i < Cell::Dim; ++i) {
    output[i] = outputEigen(i);
  }
  return output;
}

auto CellTransform::refToSpace(const std::vector<VectorT>& input) const -> std::vector<VectorT> {
  std::vector<VectorT> output(input.size());
  refToSpace(input.data(), output.data(), input.size());
  return output;
}

void CellTransform::refToSpace(const VectorT* input, VectorT* output, std::size_t count) const {
  for (std::size_t i = 0; i < count; ++i) {
    output[i] = refToSpace(input[i]);
  }
}

auto CellTransform::refToSpaceJacobianInverse(const VectorEigenT& input) const -> MatrixEigenT {
  return refToSpaceJacobian(input).inverse();
}

auto CellTransform::spaceToRefJacobian(const VectorEigenT& input) const -> MatrixEigenT {
  return refToSpaceJacobianInverse(spaceToRef(input));
}

auto CellTransform::spaceToRef(const VectorEigenT& input) const -> VectorEigenT {
  // in the general case... we need to invert a function. So... Newton.
  // we want: f(y) = x; or: f(y) - x = 0

  // start at the barycenter of the reference cell, so that the first iterate is inside the cell
  auto iterate = VectorEigenT(Cell::ReferenceBarycenter.data());

  // the residual lives in space coordinates, hence the tolerance has to scale with them
  const double eps = 1e-12 * (1.0 + input.norm());
  constexpr std::size_t Tries = 100;
  for (std::size_t i = 0; i < Tries; ++i) {
    const VectorEigenT residual = refToSpace(iterate) - input;
    if (residual.norm() < eps) {
      return iterate;
    }
    iterate -= refToSpaceJacobian(iterate).fullPivLu().solve(residual).eval();
  }

  logError() << "Root finding failed for" << input << "after" << Tries
             << "iterations. Last iterate:" << iterate << "giving" << refToSpace(iterate);
  return iterate;
}

void AffineTransform::setup(const std::array<VectorEigenT, Cell::NumVertices>& vertices) {
  offset_ = vertices[0];

  for (std::size_t i = 0; i < Cell::Dim; ++i) {
    const VectorEigenT v = vertices[i + 1] - offset_;
    for (std::size_t j = 0; j < Cell::Dim; ++j) {
      transform_(j, i) = v(j);
    }
  }

  determinant_ = transform_.determinant();
}

AffineTransform::AffineTransform(const std::array<CoordinateT, Cell::NumVertices>& vertices) {
  std::array<VectorEigenT, Cell::NumVertices> verticesEigen{};
  for (std::size_t i = 0; i < Cell::NumVertices; ++i) {
    verticesEigen[i] = VectorEigenT(vertices[i].data());
  }
  setup(verticesEigen);
}

AffineTransform::AffineTransform(const std::array<VectorEigenT, Cell::NumVertices>& vertices) {
  setup(vertices);
}

auto AffineTransform::factorization() const -> const Eigen::PartialPivLU<MatrixEigenT>& {
  if (!factorization_.has_value()) {
    factorization_.emplace(transform_);
  }
  return factorization_.value();
}

auto AffineTransform::refToSpace(const VectorEigenT& input) const -> VectorEigenT {
  return transform_ * input + offset_;
}

auto AffineTransform::refToSpaceJacobian(const VectorEigenT& /*input*/) const -> MatrixEigenT {
  // since we're linear—no dependency on the input vector here
  return transform_;
}

auto AffineTransform::fromMeshCell(std::size_t id, const MeshReader& mesh) -> AffineTransform {
  const auto& vertexIndices = mesh.getElements()[id].vertices;
  std::array<CoordinateT, Cell::NumVertices> vertices{};
  for (std::size_t i = 0; i < vertexIndices.size(); ++i) {
    vertices[i] = mesh.getVertices()[vertexIndices[i]].coords;
  }
  return AffineTransform(vertices);
}

auto AffineTransform::spaceToRef(const VectorEigenT& input) const -> VectorEigenT {
  // solving is both cheaper and more accurate than forming an explicit inverse
  return factorization().solve(input - offset_);
}

} // namespace seissol::geometry
